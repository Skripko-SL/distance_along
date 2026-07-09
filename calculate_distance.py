import argparse
import math
import csv
import time
import os

from dbfread import DBF
import numpy as np
from scipy.spatial import cKDTree

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

def load_objects(path):
    table = DBF(path, raw=True)
    field_names = list(table.field_names)
    has_id_t = 'id_t' in field_names
    objects = []
    for r in table:
        obj = {
            'id': int(r['id']),
            'lon': float(decode(r['X'])),
            'lat': float(decode(r['Y'])),
        }
        if has_id_t:
            obj['id_t'] = decode(r['id_t'])
        else:
            obj['id_t'] = str(obj['id'])
        objects.append(obj)
    return objects

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
    parser = argparse.ArgumentParser(
        description='Расчёт расстояний по прямой (great-circle) от объектов до опорной сетки'
    )
    parser.add_argument('--objects', '-o',
        default='/Users/skripko.sergey/Documents/Python/Graf/data/school.dbf',
        help='Путь к DBF с объектами (поля: id, X, Y, опционально id_t)')
    parser.add_argument('--grid', '-g',
        default='/Users/skripko.sergey/Documents/Python/Graf/data/points_buff_400.dbf',
        help='Путь к DBF с опорной сеткой')
    parser.add_argument('--output', '-O',
        default=None,
        help='Выходной CSV (по умолч.: <каталог_объектов>/<имя_объектов>_to_grid_distance.csv)')
    args = parser.parse_args()

    base_dir = os.path.dirname(args.objects)
    objects_stem = os.path.splitext(os.path.basename(args.objects))[0]
    if args.output is None:
        output_path = os.path.join(base_dir, f'{objects_stem}_to_grid_distance.csv')
    else:
        output_path = args.output

    print(f'Загрузка объектов ({objects_stem})...')
    t0 = time.time()
    objects = load_objects(args.objects)
    print(f'  Загружено {len(objects)} объектов за {time.time() - t0:.2f}с')

    print('Загрузка опорных точек сетки...')
    t0 = time.time()
    points = load_grid(args.grid)
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

    print('Поиск ближайших точек для каждого объекта...')
    t0 = time.time()
    obj_lats_rad = np.radians([s['lat'] for s in objects])
    obj_lons_rad = np.radians([s['lon'] for s in objects])
    obj_coords_3d = np.column_stack([
        np.cos(obj_lats_rad) * np.cos(obj_lons_rad),
        np.cos(obj_lats_rad) * np.sin(obj_lons_rad),
        np.sin(obj_lats_rad),
    ])
    distances_3d, indices = tree.query(obj_coords_3d, k=1)
    print(f'  Поиск выполнен за {time.time() - t0:.2f}с')

    print('Сохранение результатов...')
    t0 = time.time()

    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'source_id', 'source_id_t',
            'source_lon', 'source_lat',
            'grid_point_id', 'grid_left', 'grid_top', 'grid_col', 'grid_row',
            'distance_m'
        ])

        for i, s in enumerate(objects):
            p = points[indices[i]]
            dist_m = haversine_distance_rad(
                math.radians(s['lat']), math.radians(s['lon']),
                math.radians(p['lat']), math.radians(p['lon'])
            )

            writer.writerow([
                s['id'], s['id_t'],
                f'{s["lon"]:.8f}', f'{s["lat"]:.8f}',
                p['id'], f'{p["left"]:.2f}', f'{p["top"]:.2f}',
                p['col_index'], p['row_index'],
                f'{dist_m:.3f}'
            ])

    print(f'  Результат сохранён в {output_path} за {time.time() - t0:.2f}с')

    max_dist = max(float(r['distance_m'].replace(',', '.')) for r in csv.DictReader(open(output_path, encoding='utf-8')))
    min_dist = min(float(r['distance_m'].replace(',', '.')) for r in csv.DictReader(open(output_path, encoding='utf-8')))

    print(f'\nГотово!')
    print(f'  Всего обработано объектов: {len(objects)}')
    print(f'  Всего опорных точек: {len(points)}')
    print(f'  Мин. расстояние: {min_dist:.3f} м')
    print(f'  Макс. расстояние: {max_dist:.3f} м')

if __name__ == '__main__':
    main()
