# Инструкция по запуску расчёта расстояний и путевого времени

## Шпаргалка (быстрый старт)

**Готовые команды — просто скопировать:**

```cmd
:: школы (все файлы по умолчанию лежат рядом с .exe)
road_distance.exe

:: другой объект (больницы, поликлиники и т.п.)
road_distance.exe --objects hospitals.dbf

:: явно с учётом направленности дорог (oneway) — рекомендуемый вариант
road_distance.exe --objects school.dbf --grid all_points.dbf --roads a_graf.shp
```

**Какой файл в какой параметр класть:**

| Параметр | Что туда класть | Пример файла | Обязательные поля в DBF |
|---|---|---|---|
| `--objects` | точки объектов (школы, больницы...) | `school.dbf` | `id`, `xcoord`, `ycoord` |
| `--grid` | опорная сетка точек | `all_points.dbf` | `xcoord`, `ycoord`, `col_index`, `row_index` |
| `--roads` | дорожная сеть (линии) | `a_graf.shp` или `roads.shp` | `oneway` + геометрия `.shp` |

⚠️ **Частая ошибка:** подставить дорожный файл (`roads.dbf`/`a_graf.dbf`) в `--grid` или `--objects`, или наоборот — файл точек в `--roads`. Это разные по структуре таблицы (у дорог нет `xcoord`/`ycoord`, у точек нет `oneway`), скрипт упадёт с ошибкой чтения поля или `ShapefileException`. Всегда сверяйтесь с таблицей выше.

⚠️ **Про `--roads`:** если не указывать `--roads` явно, по умолчанию используется `a_graf.shp` — со всеми тремя файлами (`a_graf.shp`, `.dbf`, `.shx`) он должен лежать рядом с `.exe`. Если рядом лежит только `roads.*` (а `a_graf.*` нет), обязательно указывайте `--roads roads.shp`, иначе будет ошибка `Neither a_graf.dbf nor a_graf.shp could be opened`.

⚠️ **Разница `a_graf.shp` vs `roads.shp`:** обе годятся для `--roads`, но `a_graf.shp` учитывает направленность дорог (одностороннее движение), а `roads.shp` — нет (все дороги двусторонние). При сомнениях используйте `a_graf.shp`.

**Расчёт путевого времени вместо расстояния** (`road_time.exe` / `time_along_roads.py`):

```cmd
:: школы (все файлы по умолчанию лежат рядом с .exe)
road_time.exe

:: явно
road_time.exe --objects school.dbf --grid all_points.dbf --roads a_graf.shp
```

⚠️ Работает **только** с `a_graf.shp` (нужны поля `oneway` и `speed`) — `roads.shp` не подходит, скрипт/`.exe` сразу завершится с ошибкой.

---

## 1. Расчёт по дорожному графу (основной скрипт)

```bash
# школы (все пути по умолчанию)
python3 distance_along_roads.py

# больницы → достаточно указать --objects
python3 distance_along_roads.py --objects hospitals.dbf
```

Результат: `<сетка>_to_<объекты>_distance.csv` (например `all_points_to_school_distance.csv`).

### Параметры

| Флаг | Назначение | По умолчанию |
|------|-----------|-------------|
| `--objects, -o` | DBF с объектами (id, xcoord, ycoord, опц. id_t) | `school.dbf` |
| `--grid, -g` | DBF опорной сетки (xcoord, ycoord, col_index, row_index) | `all_points.dbf` |
| `--roads, -r` | Shapefile дорожной сети | `a_graf.shp` |
| `--output, -O` | Выходной CSV | `<каталог_сетки>/<сетка>_to_<объекты>_distance.csv` |
| `--k, -k` | Число кандидатов KD-дерева | 3 |

Все пути по умолчанию — в папке скрипта. Для других объектов достаточно указать `--objects`.

### Примеры

```bash
# школы (все пути по умолчанию)
python3 distance_along_roads.py

# больницы → достаточно указать --objects
python3 distance_along_roads.py --objects hospitals.dbf

# всё явно
python3 distance_along_roads.py --objects school.dbf --grid all_points.dbf --roads a_graf.shp

# свой выходной файл
python3 distance_along_roads.py -o school.dbf -O result.csv
```

### Что делает скрипт

1. Перепроецирует координаты объектов, сетки и дорог из EPSG:3857 (Web Mercator) в UTM 38N (EPSG:32638) — Web Mercator искажает расстояния в `1/cos(широта)` раз (на широте Нижегородской области — ×1.72–1.89), поэтому все геометрические расчёты ведутся в UTM, где 1 единица = 1 истинный метр
2. Строит дорожный граф из всех вершин полилиний — roads
3. Проецирует каждый объект на ближайший сегмент дороги (перпендикуляр)
4. Запускает Дейкстру от супер-источника — для каждого узла графа минимальное расстояние до любого объекта
5. Проецирует каждую точку сетки на дорогу
6. Для каждой точки: `total = перп_точки + расстояние_по_графу`
7. Сохраняет CSV

### Выходной CSV

| Колонка | Описание |
|---------|----------|
| `grid_point_id` | ID точки сетки |
| `grid_lon`, `grid_lat` | Координаты центра ячейки (WGS84) |
| `grid_col`, `grid_row` | Позиция в сетке |
| `perp_dist_m` | Перпендикуляр от точки до дороги (м) |
| `road_dist_to_source_m` | Расстояние по графу до ближайшего объекта (м) |
| `total_dist_m` | Итого: перп + граф (м) |
| `source_id` | ID ближайшего объекта |
| `source_id_t` | Строковый ID ближайшего объекта |

Недостижимые точки: `road_dist_to_source_m` = `NaN`, `total_dist_m` = `NaN`.

---

## 2. Расчёт минимального путевого времени (`time_along_roads.py`)

```bash
# школы (все пути по умолчанию)
python3 time_along_roads.py

# другой объект
python3 time_along_roads.py --objects hospitals.dbf
```

Результат: `<сетка>_to_<объекты>_time.csv` (например `all_points_to_school_time.csv`).

⚠️ В отличие от `distance_along_roads.py`, работает **только** с `a_graf.shp` — обязательны поля `oneway` и `speed` (скорость сегмента, км/ч). Если указать `--roads roads.shp` (там нет `speed`), скрипт завершится с понятной ошибкой.

### Параметры

| Флаг | Назначение | По умолчанию |
|------|-----------|-------------|
| `--objects, -o` | DBF с объектами (id, xcoord, ycoord, опц. id_t) | `school.dbf` |
| `--grid, -g` | DBF опорной сетки (xcoord, ycoord, col_index, row_index) | `all_points.dbf` |
| `--roads, -r` | Shapefile автомобильного графа (нужны `oneway`, `speed`) | `a_graf.shp` |
| `--output, -O` | Выходной CSV | `<каталог_сетки>/<сетка>_to_<объекты>_time.csv` |
| `--k, -k` | Число кандидатов KD-дерева | 3 |

### Что делает скрипт

Та же геометрия, что в `distance_along_roads.py` (перепроекция в UTM 38N, KD-дерево, перпендикуляр, брекетинг, супер-источник + Dijkstra), но веса рёбер графа — **секунды**, а не метры:

1. Перпендикуляр от объекта/точки сетки до дороги переводится в секунды по пешеходной скорости **1.4 м/с** (константа `WALK_SPEED_MPS`).
2. Движение по самой дороге (рёбра графа и частичные отрезки от узла-скобки до подошвы перпендикуляра) переводится в секунды по скорости **того сегмента дороги**, на который спроецирована точка (`speed`, км/ч → м/с).
3. `total_time_s = grid_point_time_s + network_time_s`, где `network_time_s` уже включает перпендикуляр и движение по дороге со стороны объекта.

### Выходной CSV

| Колонка | Описание |
|---------|----------|
| `grid_point_id` | ID точки сетки |
| `grid_lon`, `grid_lat` | Координаты центра ячейки (WGS84) |
| `grid_col`, `grid_row` | Позиция в сетке |
| `perp_dist_m` | Перпендикуляр от точки до дороги (м) |
| `grid_point_time_s` | Время пешком от точки до дороги (с) = `perp_dist_m / 1.4` |
| `network_time_s` | Время движения по дороге до ближайшего объекта (с) |
| `total_time_s` | Итого: `grid_point_time_s + network_time_s` (с) |
| `source_id` | ID ближайшего объекта |
| `source_id_t` | Строковый ID ближайшего объекта |

Недостижимые точки: `network_time_s` = `NaN`, `total_time_s` = `NaN`.

---

## 3. Формат входных DBF

Любой DBF-файл с полями:

| Поле | Тип | Обязательно | Описание |
|------|-----|------------|----------|
| `id` | Число | да | Уникальный ID объекта |
| `xcoord` | Число | да | X-координата (EPSG:3857, метры) |
| `ycoord` | Число | да | Y-координата (EPSG:3857, метры) |
| `id_t` | Строка | нет | Строковый идентификатор |

### Для сетки

| Поле | Тип | Описание |
|------|-----|----------|
| `xcoord` | Число | X-координата точки (EPSG:3857, метры) |
| `ycoord` | Число | Y-координата точки (EPSG:3857, метры) |
| `col_index` | Число | Индекс колонки |
| `row_index` | Число | Индекс строки |

---

## 4. Требования к среде

```bash
pip3 install scipy numpy dbfread pyshp pyproj
```

---

## 5. Сборка .exe для Windows (без Python на целевой машине)

Сборка через GitHub Actions — готовые `.exe` со всеми зависимостями внутри. Один workflow собирает сразу два файла: `road_distance.exe` (расстояния, `distance_along_roads.py`) и `road_time.exe` (путевое время, `time_along_roads.py`).

### Как получить .exe

1. Залить код на GitHub (например, `git push origin main`)
2. Перейти в репозиторий → Actions → **Build road_distance.exe и road_time.exe**
3. Скачать артефакты `road_distance` (внутри `road_distance.exe`) и `road_time` (внутри `road_time.exe`)
4. Разместить оба `.exe` и файлы данных (`*.dbf`, `a_graf.*`) в одной папке на Windows

На Windows не нужен Python и никакие библиотеки — всё упаковано в `.exe`.

### Запуск на Windows

```cmd
road_distance.exe --objects school.dbf
road_distance.exe --objects hospitals.dbf --grid all_points.dbf --roads a_graf.shp

road_time.exe --objects school.dbf
road_time.exe --objects hospitals.dbf --grid all_points.dbf --roads a_graf.shp
```

Параметры те же, что в разделах 1 и 2. `road_time.exe` требует `a_graf.shp` (поля `oneway`, `speed`) — с `roads.shp` не работает.

