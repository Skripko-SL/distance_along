import argparse
import math
import csv
import time
import os
import sys

from dbfread import DBF
import shapefile
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from pyproj import Transformer

EARTH_RADIUS_M = 6371000
MER = 20037508.34
VERTEX_PRECISION = 1

# Скорость пешехода для внедорожных отрезков (перпендикуляр от точки/объекта
# до ближайшей точки на дороге) — единая для обоих концов маршрута, м/с.
WALK_SPEED_MPS = 1.4

# Перевод скорости из км/ч (поле "speed" в a_graf) в м/с: km/h * 1000/3600.
KMH_TO_MPS = 1000.0 / 3600.0

# Исходные данные хранятся в EPSG:3857 (Web Mercator) — эта проекция сохраняет
# углы, но искажает расстояния в 1/cos(широта) раз. Все geometric-расчёты
# (перпендикуляры, рёбра графа) ведутся в UTM 38N (EPSG:32638), см.
# distance_along_roads.py / CLAUDE.md.
MERC_TO_UTM = Transformer.from_crs('EPSG:3857', 'EPSG:32638', always_xy=True)

def merc_to_utm(x, y):
    return MERC_TO_UTM.transform(x, y)

def merc_x_to_lon(x):
    return x / MER * 180.0

def merc_y_to_lat(y):
    return math.degrees(
        math.atan(math.exp(y / MER * math.pi)) * 2 - math.pi / 2
    )

def euclidean_dist(x1, y1, x2, y2):
    return math.hypot(x2 - x1, y2 - y1)

def decode(raw_bytes):
    return raw_bytes.decode('utf-8', errors='replace').strip()

def project_point_to_polyline(px, py, polyline):
    best_dist = float('inf')
    best_foot_x = px
    best_foot_y = py
    best_frac = 0.0
    total_len = 0.0
    seg_lens = []

    for i in range(len(polyline) - 1):
        x1, y1 = polyline[i]
        x2, y2 = polyline[i + 1]
        seg_len = euclidean_dist(x1, y1, x2, y2)
        seg_lens.append(seg_len)
        total_len += max(seg_len, 1e-10)

    cum_len = 0.0
    for i in range(len(polyline) - 1):
        x1, y1 = polyline[i]
        x2, y2 = polyline[i + 1]
        seg_len = seg_lens[i]
        dx, dy = x2 - x1, y2 - y1

        if seg_len < 1e-10:
            foot_x, foot_y = x1, y1
            t = 0.0
        else:
            t = ((px - x1) * dx + (py - y1) * dy) / (seg_len * seg_len)
            t = max(0.0, min(1.0, t))
            foot_x = x1 + t * dx
            foot_y = y1 + t * dy

        dist = euclidean_dist(px, py, foot_x, foot_y)
        if dist < best_dist:
            best_dist = dist
            best_foot_x = foot_x
            best_foot_y = foot_y
            best_frac = (cum_len + t * seg_len) / total_len

        cum_len += seg_len

    return best_dist, best_foot_x, best_foot_y, best_frac

def load_objects(path):
    table = DBF(path, raw=True)
    field_names = list(table.field_names)
    has_id_t = 'id_t' in field_names
    rows = list(table)
    sx_arr = np.array([float(r['xcoord']) for r in rows])
    sy_arr = np.array([float(r['ycoord']) for r in rows])
    ux_arr, uy_arr = merc_to_utm(sx_arr, sy_arr)

    objects = []
    for r, sx, sy, ux, uy in zip(rows, sx_arr, sy_arr, ux_arr, uy_arr):
        obj = {
            'id': int(r['id']),
            'lon': merc_x_to_lon(sx),
            'lat': merc_y_to_lat(sy),
            'mx': ux,
            'my': uy,
        }
        if has_id_t:
            obj['id_t'] = decode(r['id_t'])
        else:
            obj['id_t'] = str(obj['id'])
        objects.append(obj)
    return objects

def load_grid(path):
    table = DBF(path, raw=True)
    rows = list(table)
    cx_arr = np.array([float(r['xcoord']) for r in rows])
    cy_arr = np.array([float(r['ycoord']) for r in rows])
    ux_arr, uy_arr = merc_to_utm(cx_arr, cy_arr)

    points = []
    for r, cx, cy, ux, uy in zip(rows, cx_arr, cy_arr, ux_arr, uy_arr):
        points.append({
            'id': int(r['id']),
            'mx': ux,
            'my': uy,
            'lon': merc_x_to_lon(cx),
            'lat': merc_y_to_lat(cy),
            'col_index': int(r['col_index']),
            'row_index': int(r['row_index']),
        })
    return points

def build_road_time_graph(shp_path):
    """Строит граф дорог с весами рёбер во ВРЕМЕНИ (секунды), а не в метрах.

    Требует формат a_graf.shp: поле 'oneway' (направленность рёбер — 'F'/'T'/'B')
    и поле 'speed' (скорость сегмента в км/ч, одна на весь shapefile-рекорд).
    Время каждого ребра между соседними вершинами одного рекорда =
    euclidean_dist(м) / (speed(км/ч) * 1000/3600).

    Возвращает: (граф-время csr, nodes_list, segment_vert_indices, vert_tree,
    vert_to_seg, seg_speed_kmh) — последний элемент это скорость (км/ч) для
    каждого исходного shapefile-рекорда, нужна далее для перевода частичных
    отрезков "брекетинга" (от узла до подошвы перпендикуляра) в секунды.
    """
    print('  Чтение дорог...')
    t0 = time.time()
    sf = shapefile.Reader(shp_path)
    num_roads = sf.numRecords
    print(f'    Сегментов: {num_roads}')

    records = sf.records()
    field_names = [f[0] for f in sf.fields[1:]]
    if 'oneway' not in field_names or 'speed' not in field_names:
        raise SystemExit(
            f'Файл дорог "{shp_path}" не содержит полей oneway/speed — '
            f'time_along_roads.py рассчитан только на формат a_graf.shp '
            f'(автомобильный граф со скоростями)'
        )

    oneway_list = [r['oneway'] for r in records]
    seg_speed_kmh = [float(r['speed']) for r in records]

    raw_shapes = [sf.shape(i).points for i in range(num_roads)]

    print('  Перепроецирование EPSG:3857 -> UTM 38N (EPSG:32638)...')
    flat_x, flat_y, shape_offsets = [], [], [0]
    for pts in raw_shapes:
        for p in pts:
            flat_x.append(p[0])
            flat_y.append(p[1])
        shape_offsets.append(len(flat_x))
    ux_all, uy_all = merc_to_utm(np.array(flat_x), np.array(flat_y))

    coord_to_idx = {}
    nodes_list = []
    all_rows, all_cols, all_data = [], [], []
    segment_vert_indices = []
    all_vert_utm = []
    vert_to_seg = []

    total_vertices = 0
    total_edges = 0
    total_oneway_f = 0
    total_oneway_t = 0

    for i in range(num_roads):
        start, end = shape_offsets[i], shape_offsets[i + 1]
        if end - start < 2:
            segment_vert_indices.append([])
            continue

        oneway_val = oneway_list[i]
        speed_kmh = seg_speed_kmh[i]
        speed_mps = speed_kmh * KMH_TO_MPS
        if oneway_val == 'F':
            total_oneway_f += 1
        elif oneway_val == 'T':
            total_oneway_t += 1

        seg_idx_list = []
        prev_idx = None
        prev_p = None

        for j in range(start, end):
            p = (ux_all[j], uy_all[j])
            key = (round(p[0], VERTEX_PRECISION), round(p[1], VERTEX_PRECISION))
            if key not in coord_to_idx:
                coord_to_idx[key] = len(nodes_list)
                nodes_list.append(key)
            curr_idx = coord_to_idx[key]
            seg_idx_list.append(curr_idx)

            all_vert_utm.append(p)
            vert_to_seg.append(i)

            if prev_idx is not None and curr_idx != prev_idx:
                dist = euclidean_dist(p[0], p[1], prev_p[0], prev_p[1])
                if dist > 0:
                    edge_time = dist / speed_mps
                    # 'F' — только по порядку вершин (prev -> curr), 'T' —
                    # только в обратную сторону, иначе ('B') — двустороннее.
                    if oneway_val == 'F':
                        all_rows.append(prev_idx)
                        all_cols.append(curr_idx)
                        all_data.append(edge_time)
                        total_edges += 1
                    elif oneway_val == 'T':
                        all_rows.append(curr_idx)
                        all_cols.append(prev_idx)
                        all_data.append(edge_time)
                        total_edges += 1
                    else:
                        all_rows.extend([prev_idx, curr_idx])
                        all_cols.extend([curr_idx, prev_idx])
                        all_data.extend([edge_time, edge_time])
                        total_edges += 2

            prev_idx = curr_idx
            prev_p = p
            total_vertices += 1

        segment_vert_indices.append(seg_idx_list)

    num_nodes = len(nodes_list)
    print(f'    Узлов: {num_nodes}')
    print(f'    Рёбер: {total_edges}')
    print(f'    Из них односторонних сегментов: F={total_oneway_f}, T={total_oneway_t}')
    print(f'    Всего вершин (с дублями): {total_vertices}')
    print(f'    Время: {time.time()-t0:.1f}с')

    print('  Построение CSL-матрицы (веса — секунды)...')
    t0 = time.time()
    graph = csr_matrix(
        (all_data, (all_rows, all_cols)),
        shape=(num_nodes, num_nodes)
    )
    print(f'    Матрица: {num_nodes}x{num_nodes}, за {time.time()-t0:.1f}с')

    print('  Построение KD-дерева вершин...')
    t0 = time.time()
    vert_coords = np.array(all_vert_utm, dtype=np.float64)
    vert_tree = cKDTree(vert_coords)
    print(f'    {len(vert_coords)} точек, за {time.time()-t0:.1f}с')

    print('  Анализ связности...')
    t0 = time.time()
    from scipy.sparse.csgraph import connected_components
    from collections import Counter
    n_comp, labels = connected_components(csgraph=graph, directed=False)
    comp_sizes = Counter(labels)
    largest = max(comp_sizes.values())
    top5 = comp_sizes.most_common(5)
    print(f'    Компонент: {n_comp}')
    print(f'    Крупнейшая: {largest} ({largest/num_nodes*100:.1f}%)')
    for cid, sz in top5[:3]:
        print(f'      Комп {cid}: {sz} узлов')
    print(f'    Время: {time.time()-t0:.1f}с')

    return graph, nodes_list, segment_vert_indices, vert_tree, vert_to_seg, seg_speed_kmh

def find_projection(px, py, nodes_list, segment_vert_indices, vert_tree, vert_to_seg, k=3):
    dists, idxs = vert_tree.query(np.array([[px, py]]), k=k)
    checked_segs = set()
    best_dist = float('inf')
    best_info = None

    for idx in idxs[0]:
        seg_id = vert_to_seg[idx]
        if seg_id in checked_segs:
            continue
        checked_segs.add(seg_id)
        vert_indices = segment_vert_indices[seg_id]
        if len(vert_indices) < 2:
            continue

        polyline = [(nodes_list[vi][0], nodes_list[vi][1]) for vi in vert_indices]
        perp, fx, fy, frac = project_point_to_polyline(px, py, polyline)

        if perp < best_dist:
            best_dist = perp
            best_info = {
                'seg_id': seg_id,
                'perp': perp,
                'foot_x': fx,
                'foot_y': fy,
                'frac': frac,
                'vert_indices': vert_indices,
            }

    return best_info

def get_foot_bracket_vertices(frac, vert_indices, polyline):
    total_len = 0.0
    cum_lens = []
    for j in range(len(polyline) - 1):
        d = euclidean_dist(polyline[j][0], polyline[j][1], polyline[j+1][0], polyline[j+1][1])
        total_len += d
        cum_lens.append(total_len)

    if total_len <= 0:
        return vert_indices[0], vert_indices[-1], 0.0, 0.0

    target = frac * total_len
    prev_cum = 0.0
    for j, cum in enumerate(cum_lens):
        if target <= cum or j == len(cum_lens) - 1:
            n1 = vert_indices[j]
            n2 = vert_indices[j + 1]
            dist_to_n1 = max(0.0, target - prev_cum)
            dist_to_n2 = max(0.0, cum - target)
            return n1, n2, dist_to_n1, dist_to_n2
        prev_cum = cum

    n1 = vert_indices[0]
    n2 = vert_indices[-1]
    return n1, n2, 0.0, total_len

def main():
    if getattr(sys, 'frozen', False):
        script_dir = os.path.dirname(sys.executable)
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(
        description='Расчёт минимального путевого времени от объектов до опорной сетки по дорожному графу a_graf'
    )
    parser.add_argument('--objects', '-o',
        default=os.path.join(script_dir, 'school.dbf'),
        help='Путь к DBF с объектами интереса (поля: id, xcoord, ycoord, опционально id_t)')
    parser.add_argument('--grid', '-g',
        default=os.path.join(script_dir, 'all_points.dbf'),
        help='Путь к DBF с опорной сеткой')
    parser.add_argument('--roads', '-r',
        default=os.path.join(script_dir, 'a_graf.shp'),
        help='Путь к shapefile автомобильного графа (обязательны поля oneway, speed)')
    parser.add_argument('--output', '-O',
        default=None,
        help='Выходной CSV (по умолч.: <каталог_сетки>/<сетка>_to_<объекты>_time.csv)')
    parser.add_argument('--k', '-k', type=int, default=3,
        help='Количество кандидатов KD-дерева (по умолч.: 3)')
    args = parser.parse_args()

    objects_path = args.objects
    grid_path = args.grid
    roads_shp = args.roads
    K_NEAREST = args.k

    base_dir = os.path.dirname(objects_path)
    objects_stem = os.path.splitext(os.path.basename(objects_path))[0]
    if args.output is None:
        grid_stem = os.path.splitext(os.path.basename(grid_path))[0]
        output_path = os.path.join(base_dir, f'{grid_stem}_to_{objects_stem}_time.csv')
    else:
        output_path = args.output

    print('=' * 60)
    print(f'  ПУТЕВОЕ ВРЕМЯ ПО ДОРОГАМ: {objects_stem} → СЕТКА')
    print('=' * 60)
    print()

    print('[1] Загрузка данных...')
    t0 = time.time()
    objects = load_objects(objects_path)
    grid_points = load_grid(grid_path)
    print(f'  Объектов: {len(objects)}, Точек сетки: {len(grid_points)}')
    print(f'  Время: {time.time()-t0:.1f}с')
    print()

    print('[2] Построение графа времени (все вершины полилиний)...')
    graph, nodes_list, segment_vert_indices, vert_tree, vert_to_seg, seg_speed_kmh = build_road_time_graph(roads_shp)
    num_nodes = len(nodes_list)
    print()

    print(f'[3] Проекция объектов ({objects_stem}) на дороги...')
    t0 = time.time()
    source_node_min = {}
    bracket_node_source = {}
    source_proj = []
    source_id_map = {s['id']: s['id_t'] for s in objects}

    for s in objects:
        info = find_projection(s['mx'], s['my'], nodes_list, segment_vert_indices, vert_tree, vert_to_seg, K_NEAREST)
        if info is None:
            source_proj.append(None)
            continue

        vert_indices = info['vert_indices']
        polyline = [(nodes_list[vi][0], nodes_list[vi][1]) for vi in vert_indices]
        n1, n2, d1, d2 = get_foot_bracket_vertices(info['frac'], vert_indices, polyline)

        seg_speed_mps = seg_speed_kmh[info['seg_id']] * KMH_TO_MPS
        perp_time = info['perp'] / WALK_SPEED_MPS

        for n, d in [(n1, d1), (n2, d2)]:
            d_time = d / seg_speed_mps
            total_time = perp_time + d_time
            if total_time < source_node_min.get(n, float('inf')):
                source_node_min[n] = total_time
                bracket_node_source[n] = s['id']

        info['bracket_n1'] = n1
        info['bracket_n2'] = n2
        info['dist_n1'] = d1
        info['dist_n2'] = d2
        source_proj.append(info)

    print(f'  Спроецировано: {sum(1 for p in source_proj if p is not None)} из {len(objects)}')
    print(f'  Уникальных узлов доступа: {len(source_node_min)}')
    print(f'  Время: {time.time()-t0:.1f}с')
    print()

    print('[4] Супер-источник и Dijkstra (веса — секунды)...')
    t0 = time.time()

    coo = graph.tocoo()
    virt_node = num_nodes
    extended_n = num_nodes + 1
    virt_rows, virt_cols, virt_data = [], [], []
    for node, t_sec in source_node_min.items():
        virt_rows.extend([virt_node, node])
        virt_cols.extend([node, virt_node])
        virt_data.extend([t_sec, t_sec])

    all_rows = np.concatenate([coo.row, virt_rows])
    all_cols = np.concatenate([coo.col, virt_cols])
    all_data = np.concatenate([coo.data, virt_data])
    ext_graph = csr_matrix((all_data, (all_rows, all_cols)), shape=(extended_n, extended_n))

    # graph хранит рёбра в направлении, реально разрешённом движением (oneway
    # a_graf). Нужно время «от точки сетки к объекту» (против направления
    # построения дерева от супер-источника) — Dijkstra на транспонированном
    # графе с directed=True, см. distance_along_roads.py.
    dist_matrix_2d, predecessors_2d = dijkstra(
        csgraph=ext_graph.transpose().tocsr(), directed=True, indices=[virt_node],
        return_predecessors=True
    )
    node_dist = dist_matrix_2d[0, :num_nodes]
    node_pred = predecessors_2d[0, :num_nodes]
    reachable = np.isfinite(node_dist)

    print(f'  Связей от супер-источника: {len(virt_rows)}')
    print(f'  Достижимо узлов: {np.sum(reachable)} из {num_nodes} ({np.sum(reachable)/num_nodes*100:.1f}%)')
    if np.any(reachable):
        print(f'  Макс. время: {np.max(node_dist[reachable]):.0f} с')
    print(f'  Время: {time.time()-t0:.1f}с')

    print('  Распространение ID объектов по графу...')
    source_of_node = np.full(num_nodes, -1, dtype=np.int32)
    for node, obj_id in bracket_node_source.items():
        source_of_node[node] = obj_id
    sorted_nodes = np.where(reachable)[0]
    sorted_nodes = sorted_nodes[np.argsort(node_dist[sorted_nodes])]
    for node in sorted_nodes:
        pred = node_pred[node]
        if 0 <= pred < num_nodes:
            src = source_of_node[pred]
            if src >= 0:
                source_of_node[node] = src
    reachable_sources = source_of_node[reachable]
    unique_sources = len(set(reachable_sources[reachable_sources >= 0]))
    print(f'  Узлов с известным объектом: {np.sum(source_of_node >= 0)}')
    print(f'  Уникальных объектов в графе: {unique_sources}')
    print()

    print(f'[5] Проекция точек сетки...')
    t0 = time.time()
    grid_proj = []
    for gp in grid_points:
        info = find_projection(gp['mx'], gp['my'], nodes_list, segment_vert_indices, vert_tree, vert_to_seg, K_NEAREST)
        if info is not None:
            vert_indices = info['vert_indices']
            polyline = [(nodes_list[vi][0], nodes_list[vi][1]) for vi in vert_indices]
            n1, n2, d1, d2 = get_foot_bracket_vertices(info['frac'], vert_indices, polyline)
            info['bracket_n1'] = n1
            info['bracket_n2'] = n2
            info['dist_n1'] = d1
            info['dist_n2'] = d2
        grid_proj.append(info)
    print(f'  Спроецировано: {sum(1 for p in grid_proj if p is not None)} из {len(grid_points)}')
    print(f'  Время: {time.time()-t0:.1f}с')
    print()

    print('[6] Сохранение...')
    t0 = time.time()
    reachable_count = 0
    unreachable_count = 0

    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'grid_point_id', 'grid_lon', 'grid_lat',
            'grid_col', 'grid_row',
            'perp_dist_m', 'grid_point_time_s', 'network_time_s', 'total_time_s',
            'source_id', 'source_id_t',
        ])

        for i, gp in enumerate(grid_points):
            info = grid_proj[i]
            if info is None:
                writer.writerow([
                    gp['id'], f'{gp["lon"]:.8f}', f'{gp["lat"]:.8f}',
                    gp['col_index'], gp['row_index'],
                    '0.00', '0.00', 'NaN', 'NaN', '', '',
                ])
                unreachable_count += 1
                continue

            n1, n2 = info['bracket_n1'], info['bracket_n2']
            d1, d2 = info['dist_n1'], info['dist_n2']
            perp_g = info['perp']
            seg_speed_mps = seg_speed_kmh[info['seg_id']] * KMH_TO_MPS
            grid_point_time = perp_g / WALK_SPEED_MPS

            best_network_time = float('inf')
            best_sid = -1
            if reachable[n1]:
                candidate = node_dist[n1] + d1 / seg_speed_mps
                if candidate < best_network_time:
                    best_network_time = candidate
                    best_sid = source_of_node[n1]
            if reachable[n2]:
                candidate = node_dist[n2] + d2 / seg_speed_mps
                if candidate < best_network_time:
                    best_network_time = candidate
                    best_sid = source_of_node[n2]

            if best_network_time < float('inf'):
                total_time = grid_point_time + best_network_time
                reachable_count += 1
                sid = best_sid if best_sid >= 0 else ''
                sid_t = source_id_map.get(best_sid, '') if best_sid >= 0 else ''
                writer.writerow([
                    gp['id'], f'{gp["lon"]:.8f}', f'{gp["lat"]:.8f}',
                    gp['col_index'], gp['row_index'],
                    f'{perp_g:.3f}', f'{grid_point_time:.3f}',
                    f'{best_network_time:.3f}', f'{total_time:.3f}',
                    sid, sid_t,
                ])
            else:
                writer.writerow([
                    gp['id'], f'{gp["lon"]:.8f}', f'{gp["lat"]:.8f}',
                    gp['col_index'], gp['row_index'],
                    f'{perp_g:.3f}', f'{grid_point_time:.3f}', 'NaN', 'NaN', '', '',
                ])
                unreachable_count += 1

    print(f'  Сохранено в {output_path}')
    print(f'  Достижимо: {reachable_count} ({reachable_count/len(grid_points)*100:.1f}%)')
    print(f'  Недостижимо: {unreachable_count}')
    print(f'  Время: {time.time()-t0:.1f}с')
    print()
    print('ГОТОВО!')

if __name__ == '__main__':
    main()
