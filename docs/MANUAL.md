# Мануал — SynchrotronSim (командный справочник)

Все команды — из корня проекта. Интерпретатор — **`py -3.14`** (в нём numpy, pandas, pyvista, scipy,
matplotlib). Диффурное ядро — бинарь `build\SynchrotronTracker.exe`, оптимизатор —
`build\SynchrotronOptimizer.exe` (собираются из `src/` через CMake).

**Установки в БД:** `1` — СИЛА, полное кольцо (2082 магнита, замкнутое); `2` — Сибирь (сторонний
импорт, 144 магнита); `3` — SP1, первый супериод СИЛЫ (53 магнита). Эксперименты `1`/`2`/`3` ссылаются
на них же.

---

## Полный сценарий (кратко)

```powershell
# 1. Проверить .lte и, если кольцо не замкнуто, — замкнуть
py -3.14 lte_ring.py workspace/lte/configuration_v1/SYLA_desc.lte --line workspace/lte/configuration_v1/SYLA_line.lte

# 2. Загрузить конфигурацию в БД (путь к папке или файлу)
py -3.14 load_config.py workspace/lte/configuration_v1_closed --force

# 3. Трекинг
py -3.14 track3.py --mode kick-phys --exp 1 --turns 1      # матричное ядро
py -3.14 track3.py --mode diff      --exp 1                # диффурное ядро (эталон)

# 4. Посмотреть
py -3.14 visualize.py 1                                     # интерактивное окно PyVista
py -3.14 viz_html.py 1 --glyph                              # самодостаточный HTML

# 5. Зафиксировать состояние БД (она гитигнорится)
py -3.14 db_snapshot.py save --note "что это за состояние"
```

---

## 1. Проверка и замыкание .lte

**`lte_ring.py`** — универсальный инструмент, работает по ПУТИ к файлу, ничего не знает про имена
папок и установок.

```powershell
py -3.14 lte_ring.py ring.lte                                  # один файл с LINE=(...) — только отчёт
py -3.14 lte_ring.py desc.lte --line line.lte                  # последовательность отдельным файлом
py -3.14 lte_ring.py desc.lte --line line.lte --out fixed.lte  # записать замкнутую копию
py -3.14 lte_ring.py ring.lte --beamline RING --out fixed.lte  # выбрать линию, если их несколько
```

Без `--out` ничего не пишет. Форматы определяются по содержимому: голый список токенов в отдельном
файле (СИЛА), `NAME: LINE=(&...)` внутри файла (Сибирь, elegant); разворот с вложенностью, `n*ELEM`, `-SUB`.

**Что проверяет:**
- **Σ ANGLE** по кольцу — замыкание по УГЛУ (нужно ровно 2π).
- **Замыкание по ПОЛОЖЕНИЮ** — прогон всей линии **с дрейфами**, зазор конец→начало. Это ОТДЕЛЬНОЕ
  условие: угол свести мало, кольцо может не сойтись в точку.
  ⚠️ Мерить «выход последнего магнита → вход первого» **нельзя** — на стыке стоит законный длинный
  дрейф (`RING` начинается с `DR_01_2P651` = 2.651 м), и такой замер даёт метры мнимого зазора.
- **δ = E1+E2−ANGLE** — константность внутри группы (тип+длина). Если разошлась, сборка БД упадёт
  в `build_matrix_db.add_edge_delta`.

**Чем замыкает:** масштабирует ANGLE **и E1/E2** чистых диполей (без K1) на
`close_scale = (2π − Σ ANGLE_DQ) / Σ ANGLE_pure`; combined-function (K1≠0) не трогает. Края обязаны
ехать тем же множителем, иначе δ перестаёт быть константой типа.

## 2. Загрузка конфигурации в БД

**`load_config.py`** принимает ПУТЬ (папку или файл):

```powershell
py -3.14 load_config.py --list                                     # что лежит в workspace/lte/
py -3.14 load_config.py workspace/lte/configuration_v1_closed --force
py -3.14 load_config.py path/to/desc.lte --line path/to/line.lte --force
py -3.14 load_config.py path/to/ring.lte --lattice 2 --force       # сторонний плоский .lte
py -3.14 load_config.py --config v1 --force                        # короткое имя (удобство)
py -3.14 load_config.py <путь> --parse-only                        # только разбор, БД не трогается
```

| Формат | Определяется по | Пайплайн |
|---|---|---|
| **desc+line** | `*_desc.lte` + `*_line.lte` | `lte_to_lattice` → `assemble_lattice_configs` → `build_matrix_db` → `warmstart_ring` → `assemble --sync-subs` |
| **плоский** | один `.lte` с `LINE=(...)` | `s2lattice_to_db` |

- `--force` обязателен, если решётка уже есть: полная загрузка **сбрасывает `mat_*` к проектным K**
  из `.lte` (бэкапы → `workspace/tables/backups/`).
- **Конфигурации сосуществуют.** Загрузка объявляет СВОИ решётки (сама + суб-решётки из
  `sub_lattices.json`) и передаёт их вниз через `--lattices`; строки чужих установок переносятся как
  есть, со всеми столбцами. Раньше каждая загрузка СИЛЫ стирала Сибирь — теперь нет.
- **Порядок конвейера важен:** он обязан заканчиваться `assemble_lattice_configs.py --sync-subs`.
  `assemble` создаёт суб-решётки в начале, когда `mat_*` ещё пусты, а заполняют их следующие шаги и
  только у родителя — без финального шага суб-решётка остаётся с нулевой физикой.

### Отдельные шаги (вызывает load_config, можно и вручную)

| Команда | Что делает |
|---|---|
| `lte_to_lattice.py [--lte-dir DIR]` | `.lte` → `sp_elements.csv` (шаблон супериода) + `lattice_layout.csv` (геометрия экземпляров). Существующие `arg_val` сохраняются. |
| `assemble_lattice_configs.py` | Шаблон + геометрия → `lattice_configs.csv`. Запускать после правки `sp_elements.csv`. |
| `assemble_lattice_configs.py --sync-subs` | Суб-решётки = копия родителя со всеми столбцами. ПОСЛЕ warmstart. |
| `build_matrix_db.py [--lte-dir DIR] [--lattices 1,3]` | Структурные столбцы, края E1/E2 и δ из `.lte` папки `DIR`, связность. `--lattices` — какие решётки можно перезаписывать. |
| `warmstart_ring.py N [--lte-dir DIR]` | `mat_*` решётки N = проектные K из `.lte`. |
| `lte_to_design_orbit.py` | `.lte` → сегменты опорной орбиты (для оптимизатора). |
| `s2lattice_to_db.py --lte <файл> --lattice N --exp M` | Импорт стороннего плоского `.lte`. |

## 3. Откуда ядра читают конфигурацию

**Оба ядра читают `lattice_configs.csv`** — матричное (`matrix_tracker.py:186`) и диффурное C++
(`src/main.cpp:174`, по индексам столбцов). Это и есть «большая таблица»: строка на каждый экземпляр
магнита, со своей геометрией и коэффициентами. `element_connectivity.csv` задаёт **порядок** обхода.

```
.lte → lte_to_lattice → sp_elements.csv   (ШАБЛОН: слот в супериоде, правится руками)
                      → lattice_layout.csv (геометрия каждого экземпляра)
                           ↓ assemble_lattice_configs
                      → lattice_configs.csv   ← ОБА ЯДРА ЧИТАЮТ ЭТО
                           ↓ build_matrix_db  (структурные столбцы, края)
                           ↓ warmstart_ring   (mat_* = проектные K)
                           ↓ assemble --sync-subs
```

⚠️ **Адресные правки в `lattice_configs.csv` живут до следующей пересборки.** `assemble` берёт
`arg_val`/`correct_coeff` из ШАБЛОНА по ключу `(sp_type, sp_element_id)` — то есть все супериоды
получают одно значение слота, — а `mat_*` следом перезаписывает `warmstart`. Сценарий рабочий:
«собрал → поправил руками → посчитал», просто не перезагружать `.lte` после правок. Сохранение
per-instance правок между пересборками сознательно не делалось (решение 2026-07-17: не усложнять БД).

## 4. Трекинг

```powershell
py -3.14 track3.py --mode kick-phys --exp 1 [--turns N] [--thin ε_мм]
py -3.14 track3.py --mode diff      --exp 1
```

Установка и пучок берутся из `experiments.csv[exp]`; `--lattice`/`--group` переопределяют (и создают
строку опыта, если её нет). Траектории → `workspace/results/trajectories_expN.bin` — их пишет
**только `track3.py`**, поэтому после любой правки физики прогоните его заново, иначе визуализатор
покажет старое.

**Режимы:** `diff` — диффурное ядро (Boris по картам поля, эталон); `kick-phys` — матричное ядро на
физкоэф. `mat_*` из `.lte`. Токовый режим (`kick-current`) **архивирован**: связь ток↔физкоэф. не была
доверенной, ядро читает только `mat_*`.

### Параметры ядра (matrix_tracker.Settings)

Ни одного зашитого числа: приоритет **аргумент > столбец `experiments.csv` > дефолт**.

| Параметр | Столбец `experiments.csv` | Аргумент `track3.py` | Дефолт |
|---|---|---|---|
| порог потери пучка | `loss_radius_mm` | `--loss-mm` | 12 мм |
| плоскости — одна | `plane_same_tol_mm` | `--plane-same-tol-mm` | 0.1 мм |
| нормали параллельны | `plane_parallel_tol` | `--plane-parallel-tol` | 1e-9 |
| срезов тела (фолбэк) | `n_slices` | `--n-slices` | 40 |

Пустая ячейка ⇒ дефолт. `n_slices` здесь — ФОЛБЭК: приоритет у `lattice_configs[n_slices]` элемента.
Опечатка в имени параметра → `ValueError` со списком допустимых. Новые столбцы добавлять **в конец**
схемы: C++ диффур читает `experiments.csv` по индексам.

## 5. Снимок БД

`lattice_configs.csv`, `lattice_layout.csv`, `sp_elements.csv` **генерируются и гитигнорятся** — в
истории проекта их нет. Снимок закрывает эту дыру.

```powershell
py -3.14 db_snapshot.py save --note "текст" [--lte-dir <папка>]   # снять
py -3.14 db_snapshot.py check                                     # разошёлся ли с БД (код 1 — да)
py -3.14 db_snapshot.py restore --snap <папка> --force            # восстановить (с автобэкапом)
```

`check` сравнивает и профиль строк, и **содержимое**: прежний ручной снимок протух молча именно так —
размеры совпадали, а геометрия внутри была старая. Снимок → `workspace/tables/snapshot_design/` +
`manifest.json` (коммит, ветка, из какой конфигурации, профиль по решёткам).

## 6. Карты полей

| Команда | Что делает |
|---|---|
| `scan_field_maps.py` | Инвентаризация карт `workspace/fields_bin/`, `fields_raw/`. |
| `gen_synthetic_maps.py [--no-stl]` | Синтетические карты для элементов без ANSYS-карты: ресайз очищенной карты эталона под `eff_length` → `fields_bin/synthetic/`. Идемпотентно. |
| `symmetrize_maps.py` | **Драйвер:** сырой `.fld` → сглаживание + симметризация по типу мультиполя → очищенный `.bin` в `fields_bin/<тип>/`, обновляет `map_bin_path`. *(нужен scipy)* |
| `symmetrize_map.py <карта.fld> <quad\|sext\|oct\|dq> [σ]` | Одна карта вручную. |
| `clean_dipole_maps.py [dev ...]` | Медианная симметрия средней плоскости диполей (убирает паразитное медианное поле). Бэкап `.orig`. |

Карты: `.bin` (заголовок 96 б + сетка Bx,By,Bz), `.fld` (текст ANSYS). Скейл в трекере:
`map × currentScale × correct_coeff × field_direction`, `currentScale = arg_val / nominal_arg`.
Базовые карты QP/SP/OP/DQ — ANSYS 100А (`nominal_arg=100`, сетка 2 мм v3; QP_L — 5 мм v2). Диполи DP —
постоянные магниты (`nominal_arg='-'`, масштаб через `correct_coeff`).

**Число полюсов при симметризации:** quad — 4 (n=2), sext — 6 (n=3), oct — 8 (n=4). Отражения по
чётностям НЕ разделяют quad и oct — их разделяет грид-точный поворот 90°.

**Каскад при обновлении базовых карт:** новые `.fld` в `fields_raw/` + `map_raw_path`/`nominal_arg` →
`symmetrize_maps.py` → очистить `map_bin_path` у синтетики затронутых семейств + `gen_synthetic_maps.py --no-stl`
→ `clean_dipole_maps.py`. `build_matrix_db.py` при обновлении ТОЛЬКО карт не нужен (он сбрасывает `mat_*`).

## 7. Настройка

| Команда | Ядро |
|---|---|
| `optimize.py <exp> --method=global [--loss=composite\|orbit\|envelope]` | Диффур (`SynchrotronOptimizer.exe`). |
| `optimize.py <exp> --method=windowed` | Диффур, скользящее окно магнитов. |
| `optimize_combined.py <exp>` | Чередование оконного и глобального. |
| `matrix_optimizer.py warmstart` | Матричное. |
| `continuation_tuner.py --cycles N --budget B --fraction F` | Ползущий рубеж выживания по кольцу. |
| `tune.py --finish Nsp --budget B` | Одиночный проход с явным рубежом (напр. `--finish 2.5sp`). |
| `prepare_optim_config.py` | Генерация `.optimconfig` (опорный путь, окна, ручки). |

## 8. Визуализация

| Команда | Что делает |
|---|---|
| `visualize.py <exp>` | Интерактивное окно PyVista: траектория `trajectories_expN.bin` + магниты из `lattice_configs`. |
| `live_visualizer.py <exp>` | **Живой** просмотр: следит за `sp_elements.csv`, при правке — пересборка + перетрек + перерисовка. `r`/`u` — refresh, `q` — выход. |
| `viz_shot.py --exp N --mag M --view top\|side\|persp --span S --out f.png` | Off-screen скриншот с камерой на магните `M`. |
| `viz_html.py <exp> [--glyph\|--both] [--decimate 0..1] [--out f.html]` | Самодостаточный интерактивный HTML (см. §9). |
| `sync_sp_elements.py 1` | Перенести токи из `lattice_configs.csv` в `sp_elements.csv` (нужно перед `live_visualizer`). |
| `reset_magnets.py` | Сброс: все `arg_val=0`, `correct_coeff=1`. |
| `manual_kick.py` | Однократный ручной кик токов/коэф. в `sp_elements.csv`. |

**СК:** кольцо в плоскости X-Y (Z=0), вертикаль — Z, изгиб виден сверху (`--view top`). Камера
полного кольца — в `workspace/configs/visualizer.json` (орто, `parallel_scale` — полу-высота в метрах).

## 9. Передаваемый интерактивный HTML

`viz_html.py` экспортирует сцену в самодостаточный `.html` (plotly внутри файла): крутится офлайн,
передаётся файлом. Пишется в `workspace/results/html/`. Слои: магниты (цвет по типу), пучок, опорный
путь, BPM (чёрные ромбы на `kind=monitor`).

```powershell
py -3.14 viz_html.py 1 --both --out workspace/results/html/sila_ring.html
py -3.14 viz_html.py 1 --glyph            # только лёгкая (блоки) — для полного кольца
```

| Ключ | Назначение |
|---|---|
| `--both` | Обе версии: `_light` (блоки) + `_detailed` (настоящий STL). |
| `--glyph` | Только лёгкая. Для полного кольца (2082 магнита) единственный лёгкий вариант. |
| `--decimate 0..1` | Прореживание STL детальной версии. `-1` (по умолч.) — авто под ~300k треуг.; `0` форсит полный STL (для СИЛЫ — сотни МБ). |
| `--edges` / `--no-edges` | Feature-рёбра (контур формы магнитов, переключаемо в легенде). По умолч. **АВТО**: детальный полный ring — без (иначе +~27МБ), лёгкие сцены/`--glyph` — с. `--edges` форсит (СИЛА детально ⇒ ~53МБ). |
| `--out путь.html` | Своё имя (с `--both` добавятся суффиксы). |

Запускать с **`PYTHONIOENCODING=utf-8`** (иначе финальный лог падает в cp1251-консоли, HTML при этом
уже записан). Траектория авто-прореживается до ≤8000 точек/частицу.

## 10. Конвенции и грабли

- **Знак изгиба:** положительный `ANGLE` из `.lte` = поворот **ВПРАВО** (yaw убывает), как в
  elegant/MAD. Отсюда `+x` смотрит **наружу** кольца. Единственный источник конвенции —
  `lte_to_lattice.step`; `design_orbit`, `s2lattice_to_db`, `lte_ring` берут её оттуда.
  Следствие: `k1=+K1`, `k2=+K2`, `k3=+K3` из `.lte` верны дословно, а `mat_B = −ANGLE/(QP·L_eff)`
  (`warmstart_ring.py:28`). Если уложить кольцо зеркально, переворот потребуется кикам ЧЁТНЫМ по x
  (диполь ∝x⁰, секступоль ∝x²), а нечётным (квадруполь ∝x¹, октуполь ∝x³) — нет.
- **Мультиполи в изгибающих элементах** берут `u = x − x_ref`, где `x_ref` — прогиб опорной дуги
  (`Magnet._ref_bow`): ящик прямой, а опорная внутри изгиба отходит от оси на сагитту `L·θ₀/8`.
  Без этого сама опорная ловит кик `k1·сагитта` (у DQ_L1000 это ~2 мм за элемент). elegant решает то
  же конструктивно: `CSBEND` — СЕКТОР, рама гнётся вместе с опорной.
- **`sp_elements.csv` — точка правки токов руками** (шаблон супериода → всё кольцо). После правки
  нужен `assemble_lattice_configs.py`.
- **Тюн не применён намеренно:** `mat_*` = проектные K из `.lte`, чтобы сверка с elegant шла на одних
  и тех же коэффициентах.
- **Легаси** (скрипты опытов, отчёты, внешние коды) — в `legacy/`, в сценарии не нужны. Токовый режим,
  матрицы связи (`coupling_matrix`, `element_coupling`) и C++ матричное ядро (`MatrixTracker.exe`)
  архивированы — единственный матричный движок сейчас `matrix_tracker.py`.
- Сборка C++ (при изменении `src/`): CMake → `build/SynchrotronTracker.exe`, `SynchrotronOptimizer.exe`.
