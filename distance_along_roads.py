from dbfread import DBF
import shapefile
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra, connected_components
from collections import defaultdict, Counter
import math
import csv
import time
import sys

EARTH_RADIUS_M = 6371000
MER = 20037508.34
MERGE_RADIUS = 10.0

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

def build_road_graph(shp_path):
    print('  Чтение дорог...')
    t0 = time.time()
    sf = shapefile.Reader(shp_path)
    num_roads = sf.numRecords
    print(f'    Сегментов: {num_roads}')

    endpoints = []
    road_info = []
    for i in range(num_roads):
        shape = sf.shape(i)
        pts = shape.points
        if len(pts) < 2:
            continue
        p1 = (round(pts[0][0], 2), round(pts[0][1], 2))
        p2 = (round(pts[-1][0], 2), round(pts[-1][1], 2))
        rec = sf.record(i)
        leght = float(rec['leght'])
        oneway = rec['oneway']
        endpoints.append(p1)
        endpoints.append(p2)
        road_info.append((p1, p2, leght, oneway))
    print(f'    Концов: {len(endpoints)}, за {time.time()-t0:.1f}с')

    print(f'  Слияние концов в радиусе {MERGE_RADIUS}м...')
    t0 = time.time()
    eps = np.array(endpoints, dtype=np.float64)
    tree = cKDTree(eps)
    n = len(endpoints)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for idx in range(n):
        if idx % 50000 == 0:
            print(f'      {idx}/{n}', end='\r')
        neighbors = tree.query_ball_point(eps[idx], MERGE_RADIUS)
        for ni in neighbors:
            if ni > idx:
                union(idx, ni)

    print(f'      {n}/{n}')
    root_to_nodes = defaultdict(set)
    for i in range(n):
        root_to_nodes[find(i)].add(i)

    cluster_coords = {}
    for root, node_set in root_to_nodes.items():
        node_list = list(node_set)
        xs = [eps[i][0] for i in node_list]
        ys = [eps[i][1] for i in node_list]
        canonical = (round(sum(xs) / len(xs), 2), round(sum(ys) / len(ys), 2))
        cluster_coords[root] = canonical

    print(f'    Узлов после слияния: {len(cluster_coords)}, за {time.time()-t0:.1f}с')

    print('  Построение графа...')
    t0 = time.time()
    cluster_list = list(cluster_coords.values())
    cluster_to_idx = {c: i for i, c in enumerate(cluster_list)}
    cluster_tree = cKDTree(cluster_list)

    n1_list = []
    n2_list = []
    leght_list = []
    edge_set = set()

    for p1, p2, leght, oneway in road_info:
        _, i1 = cluster_tree.query(np.array([[p1[0], p1[1]]], dtype=np.float64))
        _, i2 = cluster_tree.query(np.array([[p2[0], p2[1]]], dtype=np.float64))
        n1, n2 = int(i1[0]), int(i2[0])
        if n1 == n2:
            continue
        if oneway == 'F':
            if (n1, n2) not in edge_set:
                edge_set.add((n1, n2))
                n1_list.append(n1); n2_list.append(n2); leght_list.append(leght)
        elif oneway == 'T':
            if (n2, n1) not in edge_set:
                edge_set.add((n2, n1))
                n1_list.append(n2); n2_list.append(n1); leght_list.append(leght)
        else:
            if (n1, n2) not in edge_set and (n2, n1) not in edge_set:
                edge_set.add((n1, n2))
                n1_list.append(n1); n2_list.append(n2); leght_list.append(leght)
                n1_list.append(n2); n2_list.append(n1); leght_list.append(leght)

    print(f'    Узлов: {len(cluster_list)}, Рёбер: {len(edge_set)}, за {time.time()-t0:.1f}с')

    print('  Анализ связности...')
    t0 = time.time()
    graph = csr_matrix(
        (leght_list, (n1_list, n2_list)),
        shape=(len(cluster_list), len(cluster_list))
    )
    print(f'    Матрица построена за {time.time()-t0:.1f}с')

    return graph, cluster_list, cluster_tree

def main():
    base_dir = '/Users/skripko.sergey/Documents/Python/Graf/data'
    school_path = f'{base_dir}/school.dbf'
    grid_path = f'{base_dir}/points_buff_400.dbf'
    roads_shp = f'{base_dir}/roads.shp'
    output_path = f'{base_dir}/grid_to_school_distance.csv'

    print('=' * 60)
    print('  РАСЧЁТ РАССТОЯНИЙ ПО ДОРОЖНОМУ ГРАФУ')
    print('=' * 60)
    print()

    print('[1] Загрузка данных...')
    t0 = time.time()
    schools = load_schools(school_path)
    grid_points = load_grid(grid_path)
    print(f'  Школ: {len(schools)}, Точек сетки: {len(grid_points)}, за {time.time()-t0:.1f}с')
    print()

    print('[2] Построение дорожного графа...')
    graph, cluster_list, cluster_tree = build_road_graph(roads_shp)
    num_nodes = len(cluster_list)
    print()

    print('[3] Привязка школ к графу...')
    t0 = time.time()
    school_merc = np.array([[s['mx'], s['my']] for s in schools], dtype=np.float64)
    _, school_nodes = cluster_tree.query(school_merc)
    school_nodes = school_nodes.tolist()
    print(f'  Школ привязано: {len(set(school_nodes))} уникальных узлов, за {time.time()-t0:.1f}с')
    print()

    print('[4] Мульти-источниковый Dijkstra...')
    t0 = time.time()
    sources = list(set(school_nodes))
    dist_matrix = dijkstra(
        csgraph=graph,
        directed=False,
        indices=sources,
        min_only=True,
    )
    reachable_mask = np.isfinite(dist_matrix)
    print(f'  Достижимо узлов: {np.sum(reachable_mask)} из {num_nodes} ({np.sum(reachable_mask)/num_nodes*100:.1f}%)')
    if np.any(reachable_mask):
        print(f'  Макс. расстояние: {np.max(dist_matrix[reachable_mask]):.0f} м')
    print(f'  Время: {time.time()-t0:.1f}с')
    print()

    print('[5] Привязка точек сетки к графу...')
    t0 = time.time()
    grid_merc = np.array([[gp['mx'], gp['my']] for gp in grid_points], dtype=np.float64)
    snap_dists, grid_nodes = cluster_tree.query(grid_merc)
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
            'snap_dist_m',
            'road_dist_to_school_m',
        ])

        for i, gp in enumerate(grid_points):
            node_idx = int(grid_nodes[i])
            snap_dist = float(snap_dists[i])
            d_road = float(dist_matrix[node_idx]) if reachable_mask[node_idx] else -1
            d_direct = haversine(gp['lat'], gp['lon'],
                                  schools[school_nodes.index(node_idx)]['lat'],
                                  schools[school_nodes.index(node_idx)]['lon']) if node_idx in school_nodes else haversine(
                gp['lat'], gp['lon'],
                schools[0]['lat'], schools[0]['lon']
            )

            if d_road >= 0:
                reachable_count += 1
            else:
                unreachable_count += 1

            writer.writerow([
                gp['id'],
                f'{gp["lon"]:.8f}',
                f'{gp["lat"]:.8f}',
                gp['col_index'],
                gp['row_index'],
                f'{snap_dist:.2f}',
                f'{d_road:.3f}' if d_road >= 0 else 'NaN',
            ])

    print(f'  Сохранено в {output_path} за {time.time()-t0:.1f}с')
    print(f'  Достижимо: {reachable_count} ({reachable_count/len(grid_points)*100:.1f}%)')
    print(f'  Недостижимо: {unreachable_count} ({unreachable_count/len(grid_points)*100:.1f}%)')
    print()
    print('ГОТОВО!')

if __name__ == '__main__':
    main()
