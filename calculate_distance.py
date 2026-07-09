from dbfread import DBF
import numpy as np
from scipy.spatial import cKDTree
import math
import csv
import sys
import time

EARTH_RADIUS_M = 6371000

def merc_x_to_lon(x):
    return x / 20037508.34 * 180.0

def merc_y_to_lat(y):
    return math.degrees(
        math.atan(math.exp(y / 20037508.34 * math.pi)) * 2 - math.pi / 2
    )

def haversine_distance_rad(lat1_rad, lon1_rad, lat2_rad, lon2_rad):
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return c * EARTH_RADIUS_M

def decode(raw_bytes):
    return raw_bytes.decode('utf-8', errors='replace').strip()

def load_schools(path):
    table = DBF(path, raw=True)
    schools = []
    for r in table:
        schools.append({
            'id': int(r['id']),
            'id_t': decode(r['id_t']),
            'lon': float(decode(r['X'])),
            'lat': float(decode(r['Y'])),
            'name': decode(r['name']),
            'addres': decode(r['addres']),
        })
    return schools

def load_grid(path):
    table = DBF(path, raw=True)
    points = []
    for r in table:
        left = float(r['left'])
        top = float(r['top'])
        cx = left + 200
        cy = top + 200
        lon = merc_x_to_lon(cx)
        lat = merc_y_to_lat(cy)
        points.append({
            'id': int(r['id']),
            'left': left,
            'top': top,
            'col_index': int(r['col_index']),
            'row_index': int(r['row_index']),
            'lon': lon,
            'lat': lat,
        })
    return points

def main():
    base_dir = '/Users/skripko.sergey/Documents/Python/Graf/data'
    school_path = f'{base_dir}/school.dbf'
    grid_path = f'{base_dir}/points_buff_400.dbf'
    output_path = f'{base_dir}/school_to_grid_distance.csv'

    print('Загрузка школ...')
    t0 = time.time()
    schools = load_schools(school_path)
    print(f'  Загружено {len(schools)} школ за {time.time() - t0:.2f}с')

    print('Загрузка опорных точек сетки...')
    t0 = time.time()
    points = load_grid(grid_path)
    print(f'  Загружено {len(points)} точек за {time.time() - t0:.2f}с')

    print('Построение KD-дерева (3D-декартовы координаты)...')
    t0 = time.time()
    lats_rad = np.radians([p['lat'] for p in points])
    lons_rad = np.radians([p['lon'] for p in points])
    coords_3d = np.column_stack([
        np.cos(lats_rad) * np.cos(lons_rad),
        np.cos(lats_rad) * np.sin(lons_rad),
        np.sin(lats_rad),
    ])
    tree = cKDTree(coords_3d)
    print(f'  Дерево построено за {time.time() - t0:.2f}с')

    print('Поиск ближайших точек для каждой школы...')
    t0 = time.time()
    school_lats_rad = np.radians([s['lat'] for s in schools])
    school_lons_rad = np.radians([s['lon'] for s in schools])
    school_coords_3d = np.column_stack([
        np.cos(school_lats_rad) * np.cos(school_lons_rad),
        np.cos(school_lats_rad) * np.sin(school_lons_rad),
        np.sin(school_lats_rad),
    ])
    distances_3d, indices = tree.query(school_coords_3d, k=1)
    print(f'  Поиск выполнен за {time.time() - t0:.2f}с')

    print('Сохранение результатов...')
    t0 = time.time()

    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'school_id', 'school_id_t', 'school_name', 'school_address',
            'school_x', 'school_y',
            'grid_point_id', 'grid_left', 'grid_top', 'grid_col', 'grid_row',
            'distance_m'
        ])

        for i, s in enumerate(schools):
            p = points[indices[i]]
            dist_m = haversine_distance_rad(
                math.radians(s['lat']), math.radians(s['lon']),
                math.radians(p['lat']), math.radians(p['lon'])
            )

            writer.writerow([
                s['id'], s['id_t'], s['name'], s['addres'],
                f'{s["lon"]:.8f}', f'{s["lat"]:.8f}',
                p['id'], f'{p["left"]:.2f}', f'{p["top"]:.2f}',
                p['col_index'], p['row_index'],
                f'{dist_m:.3f}'
            ])

    print(f'  Результат сохранён в {output_path} за {time.time() - t0:.2f}с')

    max_dist = max(float(r['distance_m'].replace(',', '.')) for r in csv.DictReader(open(output_path, encoding='utf-8')))
    min_dist = min(float(r['distance_m'].replace(',', '.')) for r in csv.DictReader(open(output_path, encoding='utf-8')))

    print(f'\nГотово!')
    print(f'  Всего обработано школ: {len(schools)}')
    print(f'  Всего опорных точек: {len(points)}')
    print(f'  Мин. расстояние до ближайшей точки: {min_dist:.3f} м')
    print(f'  Макс. расстояние до ближайшей точки: {max_dist:.3f} м')

if __name__ == '__main__':
    main()
