from dbfread import DBF
import shapefile
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix, vstack, hstack
from scipy.sparse.csgraph import dijkstra
from collections import defaultdict
import math
import csv
import time

EARTH_RADIUS_M = 6371000
MER = 20037508.34
NODE_PRECISION = 0

def merc_x_to_lon(x):
    return x / MER * 180.0

def merc_y_to_lat(y):
    return math.degrees(
        math.atan(math.exp(y / MER * math.pi)) * 2 - math.pi / 2
    )

def haversine(lat1, lon1, lat2, lon2):
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat / 2) ** 2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_M * math.atan2(math.sqrt(a), math.sqrt(1 - a))

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
        dx, dy = x2 - x1, y2 - y1
        seg_len = math.hypot(dx, dy)
        seg_lens.append(seg_len)
        total_len += seg_len if seg_len > 0 else 1e-10

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

        dist = math.hypot(px - foot_x, py - foot_y)
        if dist < best_dist:
            best_dist = dist
            best_foot_x = foot_x
            best_foot_y = foot_y
            seg_frac = (cum_len + t * seg_len) / total_len if total_len > 0 else 0.0
            best_frac = seg_frac

        cum_len += seg_len

    return best_dist, best_foot_x, best_foot_y, best_frac

def load_schools(path):
    table = DBF(path, raw=True)
    schools = []
    for r in table:
        sx = float(r['X']) / 180 * MER
        sy = MER / math.pi * math.log(math.tan(math.pi / 4 + math.radians(float(r['Y'])) / 2))
        schools.append({
            'id': int(r['id']),
            'id_t': decode(r['id_t']),
            'lon': float(decode(r['X'])),
            'lat': float(decode(r['Y'])),
            'name': decode(r['name']),
            'addres': decode(r['addres']),
            'mx': sx,
            'my': sy,
        })
    return schools

def load_grid(path):
    table = DBF(path, raw=True)
    points = []
    for r in table:
        cx = float(r['left']) + 200
        cy = float(r['top']) + 200
        lon = merc_x_to_lon(cx)
        lat = merc_y_to_lat(cy)
        points.append({
            'id': int(r['id']),
            'mx': cx,
            'my': cy,
            'lon': lon,
            'lat': lat,
            'col_index': int(r['col_index']),
            'row_index': int(r['row_index']),
        })
    return points

def build_road_graph_and_index(shp_path):
    print('  Чтение дорог...')
    t0 = time.time()
    sf = shapefile.Reader(shp_path)
    num_roads = sf.numRecords
    print(f'    Сегментов: {num_roads}')

    coord_to_idx = {}
    nodes_list = []
    edges = []
    segments = []
    all_vert_coords = []
    vert_to_seg = []

    for i in range(num_roads):
        shape = sf.shape(i)
        rec = sf.record(i)
        pts = shape.points
        if len(pts) < 2:
            continue

        leght = float(rec['leght'])
        oneway = rec['oneway']

        k1 = (round(pts[0][0], NODE_PRECISION), round(pts[0][1], NODE_PRECISION))
        k2 = (round(pts[-1][0], NODE_PRECISION), round(pts[-1][1], NODE_PRECISION))

        for k in (k1, k2):
            if k not in coord_to_idx:
                coord_to_idx[k] = len(nodes_list)
                nodes_list.append(k)

        n1, n2 = coord_to_idx[k1], coord_to_idx[k2]

        if oneway == 'F':
            edges.append((n1, n2, leght))
        elif oneway == 'T':
            edges.append((n2, n1, leght))
        else:
            edges.append((n1, n2, leght))
            edges.append((n2, n1, leght))

        seg_pts = [(p[0], p[1]) for p in pts]
        seg_id = len(segments)
        segments.append({
            'id': seg_id,
            'points': seg_pts,
            'n1': n1,
            'n2': n2,
            'length': leght,
            'oneway': oneway,
        })

        for p in pts:
            all_vert_coords.append((p[0], p[1]))
            vert_to_seg.append(seg_id)

    num_nodes = len(nodes_list)
    print(f'    Узлов: {num_nodes}, сегментов: {len(segments)}')
    print(f'    Всего вершин полилиний: {len(all_vert_coords)}')
    print(f'    Время: {time.time()-t0:.1f}с')

    print('  Построение CSL-матрицы...')
    t0 = time.time()
    row = [e[0] for e in edges]
    col = [e[1] for e in edges]
    data = [e[2] for e in edges]
    graph = csr_matrix((data, (row, col)), shape=(num_nodes, num_nodes))
    print(f'    Матрица: {num_nodes}x{num_nodes}, рёбер: {len(edges)}, за {time.time()-t0:.1f}с')

    print('  Построение KD-дерева вершин...')
    t0 = time.time()
    vert_coords = np.array(all_vert_coords, dtype=np.float64)
    vert_tree = cKDTree(vert_coords)
    print(f'    {len(vert_coords)} точек, за {time.time()-t0:.1f}с')

    return graph, nodes_list, segments, vert_tree, vert_to_seg

def find_projection(px, py, segments, vert_tree, vert_to_seg, k=3):
    dists, idxs = vert_tree.query(np.array([[px, py]]), k=k)
    checked_segs = set()
    best_dist = float('inf')
    best_info = None

    for idx in idxs[0]:
        seg_id = vert_to_seg[idx]
        if seg_id in checked_segs:
            continue
        checked_segs.add(seg_id)
        seg = segments[seg_id]
        perp, fx, fy, frac = project_point_to_polyline(px, py, seg['points'])

        if perp < best_dist:
            best_dist = perp
            best_info = {
                'seg_id': seg_id,
                'perp': perp,
                'foot_x': fx,
                'foot_y': fy,
                'frac': frac,
                'n1': seg['n1'],
                'n2': seg['n2'],
                'length': seg['length'],
            }

    return best_info

def main():
    base_dir = '/Users/skripko.sergey/Documents/Python/Graf/data'
    school_path = f'{base_dir}/school.dbf'
    grid_path = f'{base_dir}/points_buff_400.dbf'
    roads_shp = f'{base_dir}/roads.shp'
    output_path = f'{base_dir}/grid_to_school_distance.csv'
    K_NEAREST = 3

    print('=' * 60)
    print('  РАССТОЯНИЕ ПО ДОРОГАМ: ПЕРПЕНДИКУЛЯР + ГРАФ + ПЕРПЕНДИКУЛЯР')
    print('=' * 60)
    print()

    print('[1] Загрузка данных...')
    t0 = time.time()
    schools = load_schools(school_path)
    grid_points = load_grid(grid_path)
    print(f'  Школ: {len(schools)}, Точек сетки: {len(grid_points)}')
    print(f'  Время: {time.time()-t0:.1f}с')
    print()

    print('[2] Построение графа и индекса...')
    graph, nodes_list, segments, vert_tree, vert_to_seg = build_road_graph_and_index(roads_shp)
    num_nodes = len(nodes_list)
    print()

    print(f'[3] Проекция школ на дороги (перпендикуляр)...')
    t0 = time.time()
    school_node_min = {}  # node -> min distance from school to this node via any school
    school_proj = []

    for s in schools:
        info = find_projection(s['mx'], s['my'], segments, vert_tree, vert_to_seg, K_NEAREST)
        if info is None:
            school_proj.append(None)
            continue

        frac = info['frac']
        L = info['length']
        dist_to_n1 = frac * L
        dist_to_n2 = (1 - frac) * L

        school_node_min[info['n1']] = min(
            school_node_min.get(info['n1'], float('inf')),
            info['perp'] + dist_to_n1
        )
        school_node_min[info['n2']] = min(
            school_node_min.get(info['n2'], float('inf')),
            info['perp'] + dist_to_n2
        )

        school_proj.append(info)

    n_schools_projected = sum(1 for p in school_proj if p is not None)
    n_school_nodes = len(school_node_min)
    print(f'  Спроецировано школ: {n_schools_projected} из {len(schools)}')
    print(f'  Уникальных узлов с доступом к школе: {n_school_nodes}')
    print(f'  Время: {time.time()-t0:.1f}с')
    print()

    print('[4] Построение супер-источника и запуск Dijkstra...')
    t0 = time.time()

    # Создаём расширенную матрицу с виртуальным супер-узлом
    virt_node = num_nodes
    extended_n = num_nodes + 1

    coo = graph.tocoo()
    virt_rows, virt_cols, virt_data = [], [], []

    for node, dist in school_node_min.items():
        virt_rows.append(virt_node)
        virt_cols.append(node)
        virt_data.append(dist)
        virt_rows.append(node)
        virt_cols.append(virt_node)
        virt_data.append(dist)

    if virt_rows:
        all_rows = np.concatenate([coo.row, virt_rows])
        all_cols = np.concatenate([coo.col, virt_cols])
        all_data = np.concatenate([coo.data, virt_data])
        extended_graph = csr_matrix(
            (all_data, (all_rows, all_cols)),
            shape=(extended_n, extended_n)
        )
    else:
        extended_graph = csr_matrix((extended_n, extended_n), dtype=np.float64)

    dist_matrix = dijkstra(
        csgraph=extended_graph,
        directed=False,
        indices=[virt_node],
        min_only=True,
    )

    node_dist = dist_matrix[:num_nodes]
    reachable = np.isfinite(node_dist)
    print(f'  Супер-источник: {len(virt_rows)} связей')
    print(f'  Достижимо узлов: {np.sum(reachable)} из {num_nodes} ({np.sum(reachable)/num_nodes*100:.1f}%)')
    if np.any(reachable):
        print(f'  Макс. расстояние: {np.max(node_dist[reachable]):.0f} м')
    print(f'  Время: {time.time()-t0:.1f}с')
    print()

    print(f'[5] Проекция точек сетки на дороги (перпендикуляр)...')
    t0 = time.time()
    grid_proj = []
    for gp in grid_points:
        info = find_projection(gp['mx'], gp['my'], segments, vert_tree, vert_to_seg, K_NEAREST)
        grid_proj.append(info)
    n_grid_projected = sum(1 for p in grid_proj if p is not None)
    print(f'  Спроецировано точек: {n_grid_projected} из {len(grid_points)}')
    print(f'  Время: {time.time()-t0:.1f}с')
    print()

    print('[6] Сохранение результатов...')
    t0 = time.time()
    reachable_count = 0
    unreachable_count = 0

    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'grid_point_id', 'grid_lon', 'grid_lat',
            'grid_col', 'grid_row',
            'perp_dist_m',
            'road_dist_to_school_m',
            'total_dist_m',
        ])

        for i, gp in enumerate(grid_points):
            g_info = grid_proj[i]

            if g_info is None:
                writer.writerow([
                    gp['id'],
                    f'{gp["lon"]:.8f}', f'{gp["lat"]:.8f}',
                    gp['col_index'], gp['row_index'],
                    '0.00', 'NaN', 'NaN',
                ])
                unreachable_count += 1
                continue

            perp_g = g_info['perp']
            frac_g = g_info['frac']
            L_g = g_info['length']
            n1_g, n2_g = g_info['n1'], g_info['n2']

            d_n1_g = frac_g * L_g
            d_n2_g = (1 - frac_g) * L_g

            road_dist = float('inf')
            best_school_perp = 0

            if reachable[n1_g]:
                d = node_dist[n1_g] + d_n1_g
                if d < road_dist:
                    road_dist = d
            if reachable[n2_g]:
                d = node_dist[n2_g] + d_n2_g
                if d < road_dist:
                    road_dist = d

            if road_dist < float('inf'):
                total = perp_g + road_dist
                reachable_count += 1
                writer.writerow([
                    gp['id'],
                    f'{gp["lon"]:.8f}', f'{gp["lat"]:.8f}',
                    gp['col_index'], gp['row_index'],
                    f'{perp_g:.3f}',
                    f'{road_dist:.3f}',
                    f'{total:.3f}',
                ])
            else:
                writer.writerow([
                    gp['id'],
                    f'{gp["lon"]:.8f}', f'{gp["lat"]:.8f}',
                    gp['col_index'], gp['row_index'],
                    f'{perp_g:.3f}',
                    'NaN', 'NaN',
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
