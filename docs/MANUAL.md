# Мануал — SynchrotronSim (командный справочник)

Все команды — из корня проекта. Интерпретатор — **`py -3.14`** (в нём есть numpy, pandas,
pyvista, scipy, matplotlib). Альтернатива для трекинга без визуализации/scipy —
`.\venv\Scripts\python.exe`. Диффурное ядро — бинарь `build\SynchrotronTracker.exe`,
оптимизатор — `build\SynchrotronOptimizer.exe` (собираются из `src/` через CMake).

Значения по умолчанию ниже: решётка **101** (полное кольцо), группа пучка **1**.

---

## Полный сценарий (кратко)

```powershell
# 1. LTE → БД
py -3.14 lte_to_lattice.py                 # .lte → sp_elements.csv + lattice_layout.csv
py -3.14 assemble_lattice_configs.py       # → lattice_configs.csv (полное кольцо)
py -3.14 lte_warmstart.py                  # warm-start физкоэф. mat_* из проектных K
py -3.14 build_matrix_db.py                # доп. таблицы матричной модели

# 2. Матричный трекинг (на физкоэф.) — проверить замыкание кольца
py -3.14 track3.py --mode kick-phys --lattice 101 --group 1 --exp 60

# 3. Матрицы связи (сравнение ядер) → переход к токам
py -3.14 build_coupling_pipeline.py        # собрать coupling_matrix.csv
py -3.14 phys_to_current.py 101            # mat_* → токи (arg_val/correct_coeff)

# 4. Диффурный трекинг (на токах) — подтвердить расчёт
py -3.14 track3.py --mode diff --lattice 101 --group 1 --exp 60

# 5. Траектории и BPM выгружаются в БД автоматически (см. §6)

# 6. Настройка токов (если пучок теряется)
py -3.14 continuation_tuner.py --cycles 5 --budget 5000 --fraction 1.0
```

---

## 1. Загрузка LTE и заполнение БД

### Единый загрузчик конфигураций (рекомендуется)

Конфигурации решёток лежат в подпапках `workspace/lte/` (`configuration_v1`, `configuration_v3_syberia`, …).
**`load_config.py`** грузит любую из них ОДНОЙ командой, читая `.lte` ПРЯМО из папки (без копий в корень) и
авто-определяя формат:

```powershell
py -3.14 load_config.py --list                    # показать конфигурации и их формат
py -3.14 load_config.py --config v1               # СИЛА: полный пайплайн → готовая решётка
py -3.14 load_config.py --config syberia          # сторонний .lte → s2lattice (реш. 200)
py -3.14 load_config.py --config v1 --parse-only  # только разбор в таблицы (недеструктивно)
```

| Формат папки | Определяется по | Пайплайн (folder-aware) |
|---|---|---|
| **СИЛА** | `*_desc.lte` + `*_line.lte` | `lte_to_lattice` → `assemble` → `build_matrix_db` → `warmstart_ring` |
| **плоский** (PyAT/elegant) | один `.lte` с `LINE=(...)` | `s2lattice_to_db` (геометрия + `mat_*`) |
| **none** | нет `.lte` (только xlsx/pdf) | пропуск |

- `--config <имя>` — подстрока имени папки (`v1`, `syberia`, `v3`).
- `--force` — перезаписать решётку, которая УЖЕ есть в БД (иначе загрузчик останавливается — защита).
  Полная загрузка СИЛА **сбрасывает тюн**: `build_matrix_db` перезапишет откалиброванную `coupling_matrix`,
  `warmstart` — `mat_*` проектными K (бэкапы → `workspace/tables/backups/`). После — `recalibrate_coupling` → `phys_to_current`.
- `--parse-only` — остановиться после разбора (`sp_elements`+`lattice_layout` для СИЛА), не трогая `lattice_configs`/тюн.
- `--lattice N` — переопределить id решётки; `--gen-assets` — генерировать синтетические карты/STL (СИЛА).

### Отдельные шаги (вызываются load_config, можно и вручную)

| Команда | Что делает |
|---|---|
| `py -3.14 lte_to_lattice.py [--lte-dir DIR]` | Парсит `*_desc.lte`+`*_line.lte` (из `DIR`, по умолч. `workspace/lte`) → `sp_elements.csv` (шаблон супериода) + `lattice_layout.csv` (развёртка кольца). Существующие `arg_val` сохраняются. |
| `py -3.14 assemble_lattice_configs.py` | Разворачивает шаблон по кольцу → `lattice_configs.csv`. **Запускать после любой правки `sp_elements.csv`.** |
| `py -3.14 lte_warmstart.py` | Warm-start физкоэф. `mat_*` из проектных `K` (.lte). |
| `py -3.14 warmstart_ring.py [N] [--lte-dir DIR]` | Warm-start всего кольца N (standard + junction) проектными K из `.lte` в `DIR`. |
| `py -3.14 build_matrix_db.py` | Сборка вспомогательных таблиц матричной модели. |
| `py -3.14 lte_to_design_orbit.py` | .lte → сегменты опорной орбиты (для оптимизатора). |

## 2. Карты полей

| Команда | Что делает |
|---|---|
| `py -3.14 scan_field_maps.py` | Инвентаризация карт `workspace/fields_bin/`, `fields_raw/`. |
| `py -3.14 gen_synthetic_maps.py [--no-stl]` | Генерирует синтетические карты для элементов библиотеки без ANSYS-карты, ресайзом ОЧИЩЕННОЙ карты эталона (`map_bin_path`) под `eff_length` → `fields_bin/synthetic/`. Идемпотентно (пропускает элементы с картой). |

Карты: `.bin` (заголовок 96 б + сетка Bx,By,Bz), `.fld` (текст ANSYS). Базовые карты разложены по
подпапкам `fields_bin/<тип>/` (`quadrupole/ sextupole/ octupole/ dipole_quadrupole/ dipole/ synthetic/`);
ссылка — `hardware_library.map_bin_path` (относит. `fields_bin/`). Скейл поля в трекере:
`map × currentScale × correct_coeff × field_direction`, `currentScale = arg_val / nominal_arg`.
**Базовые карты QP/SP/OP/DQ — пересчитанные ANSYS 100А (`nominal_arg=100`, 2 мм сетка v3;
QP_L — 5 мм v2)**, лежат в `fields_raw/` (`*_100A_*.fld`). Поле 100А-абсолют → `nominal_arg=100`
и у базовых, и у синтетики (иначе ×100). Диполи DP — постоянные магниты (`nominal_arg='-'`, масштаб
через `correct_coeff`). *Прежние `*_1A` карты (утечка ток→сердечник) — в архиве, из БД убраны.*

## 3. Сглаживание и симметризация карт

Симметризация = проекция карты на симметрию мультиполя (убирает шум/асимметрию неверной симметрии,
сохраняя доминантный мультиполь). **Важно учитывать число полюсов**: quad — 4 полюса (n=2), sext — 6
(n=3), oct — 8 (n=4). Отражения по чётностям НЕ разделяют quad и oct (одинаковые чётности) — их
разделяет грид-точный поворот 90° (на квадратной сетке).

| Команда | Что делает |
|---|---|
| `py -3.14 symmetrize_maps.py` | **Драйвер (рекомендуется):** читает сырой `.fld` каждого базового элемента из `fields_raw/` (по `map_raw_path`), сглаживает+симметризует с учётом типа: quad/sext/oct → правильная мультипольная симметрия; DQ (combined-function) → только сглаживание; DP → перенос. Пишет очищенный `.bin` с `nominalArg` (из `nominal_arg`) в подпапки `fields_bin/<тип>/`, обновляет `hardware_library.map_bin_path` (бэкап `.bak`). *(нужен scipy)* |
| `py -3.14 symmetrize_map.py <карта.fld> <quad\|sext\|oct\|dq> [σ]` | Одна карта: проекция на симметрию мультиполя (отражения + поворот 90° для quad/oct) + Гаусс-сглаживание σ узлов с сохранением энергии. Пишет `_sym.bin` и `_symsm.bin`. |
| `py -3.14 clean_dipole_maps.py [dev ...]` | Медианная симметрия средней плоскости (Bx,By нечёт / Bz чёт по z), убрать паразитное медианное поле диполей. Без аргументов — все `DP_*` в `synthetic/`. Бэкап `.orig`. |

**Порядок при обновлении базовых карт** (полный каскад): (1) скопировать новые `.fld` в `fields_raw/`,
прописать `map_raw_path` + `nominal_arg` (для 100А-карт =100, у базовых И синтетики семейства) →
(2) `symmetrize_maps.py` (очищенные `.bin` из `.fld`) → (3) очистить `map_bin_path` у синтетики
изменённых семейств и `gen_synthetic_maps.py --no-stl` (пересоздать из очищенных эталонов) →
(4) `clean_dipole_maps.py` (DP-синтетика) → (5) `recalibrate_coupling.py` (матрица связи вслед за
картами; квадруполи пробуются при рабочем ~1А — сильные новые карты при 10А дают толстую линзу и
смещают `per_G`). `build_matrix_db.py` при обновлении ТОЛЬКО карт **не нужен** (он сбрасывает тюн `mat_*`).

## 4. Матричное (киковое) ядро

| Команда | Что делает |
|---|---|
| `py -3.14 track3.py --mode kick-phys --lattice 101 --group 1 --exp N` | Трекинг на **проектных физкоэф.** `mat_*`. Пишет `workspace/results/trajectories_expN.bin`. |
| `py -3.14 matrix_traj_export.py [exp] [lattice]` | Экспорт траектории матричной модели в формат визуализатора. |

## 5. Сравнение моделей и матрицы связи

| Команда | Что делает |
|---|---|
| `py -3.14 build_coupling_pipeline.py` | Собирает `coupling_matrix.csv` (ток→физкоэф.) из достоверных источников. `--dry-run` — только отчёт. |
| `py -3.14 combo_coupling.py` | Калибровка матриц связи многообмоточных (SP, DQ) по КОМБИНАЦИЯМ токов. |
| `py -3.14 element_matrix_compare.py` | Сравнение траекторий одного элемента: Boris по картам vs матричная модель. |
| `py -3.14 validate_coupling_table.py` | Валидация `element_coupling.csv` через матричный трекер. |
| `py -3.14 single_element_boris.py` | Эталонный прогон ОДНОГО элемента через Boris. |
| `py -3.14 phys_to_current.py 101` | Обратная задача: `mat_*` → токи (для запуска на токах). Флаг `--dipole-priority` — приоритет изгиба у вырожденного DQ (иначе min-norm). Бэкап `lattice_configs.csv.bak`. |

## 6. Диффурное ядро и 3 режима трекинга

| Команда | Режим |
|---|---|
| `py -3.14 track3.py --mode diff --lattice 101 --group 1 --exp N` | **Диффур** (Boris по картам, на токах). Эталон. Запускает `SynchrotronTracker.exe`. |
| `py -3.14 track3.py --mode kick-phys ...` | Матрица на физкоэф. `mat_*`. |
| `py -3.14 track3.py --mode kick-current ...` | Матрица на токах → `coupling_matrix` → `mat_*`. |

**Выгрузка в БД:**
- **Траектории** → `workspace/results/trajectories_expN.bin` (пишется автоматически трекером).
- **Показания BPM** → `workspace/tables/bpm_readings.csv`. Пишутся C++ ядром автоматически,
  если в решётке есть элементы `kind=monitor` и задан `bpm_config.csv` (несовершенства
  датчиков: dx/dy/tilt/gain/setpoint). Оба ядра (C++ и Python) дают идентичные показания.

## 7. Настройка токов

### Глобальная (все ручки сразу)

| Команда | Ядро |
|---|---|
| `py -3.14 optimize.py <exp> --method=global [--loss=composite\|orbit\|envelope]` | Диффур (`SynchrotronOptimizer.exe`). |
| `py -3.14 matrix_optimizer.py warmstart` | Матричное. |
| `py -3.14 continuation_tuner.py --cycles N --budget B --fraction F` | Авто-продолжатель: ползущий рубеж выживания по кольцу (диффур). |
| `py -3.14 tune.py --finish Nsp --budget B` | Одиночный проход с явным рубежом (напр. `--finish 2.5sp`). |

### Оконная (скользящее окно магнитов)

| Команда | Что делает |
|---|---|
| `py -3.14 optimize.py <exp> --method=windowed` | Оконная оптимизация на диффуре. |
| `py -3.14 optimize_combined.py <exp>` | Чередование оконного и глобального. |
| `py -3.14 prepare_optim_config.py` | Генерация `.optimconfig` (опорный путь, окна, ручки). |
| `py -3.14 make_nsp_drift_configs.py` / `make_nsp_next_configs.py` | Конфиги «N супериодов + дрейф / начало SP(N+1)». |

## 8. Визуализация и живая настройка

| Команда | Что делает |
|---|---|
| `py -3.14 live_visualizer.py <exp>` | **Живой** просмотр: окно PyVista следит за `sp_elements.csv`; сохранил правку токов → авто-пересборка + перетрек диффуром + перерисовка. Клавиши `r`/`u` — ручной refresh, `q` — выход. |
| `py -3.14 visualize.py <exp>` | Статичный просмотр траектории `trajectories_expN.bin`. |
| `py -3.14 viz_shot.py --exp N --mag M --view top\|side\|persp --span S --out f.png` | Off-screen скриншот сцены с камерой на магните `M`. |
| `py -3.14 viz_html.py <exp> [--glyph] [--decimate 0..1] [--out f.html]` | **Экспорт сцены в самодостаточный интерактивный HTML** (магниты+пучок+опорный путь+BPM). Крутится/зумится в любом браузере офлайн, передаётся файлом. → `workspace/results/html/`. См. §9. |
| `py -3.14 sync_sp_elements.py 101` | Перенести настроенные токи из `lattice_configs.csv` в `sp_elements.csv` (нужно перед `live_visualizer`, т.к. тюн пишет в `lattice_configs`). |
| `py -3.14 reset_magnets.py` | Сброс к дефолту: все `arg_val=0`, `correct_coeff=1`. |
| `py -3.14 manual_kick.py` | Однократный ручной кик токов/коэф. в `sp_elements.csv`. |

**Обзор всего кольца:** камера настроена в `workspace/configs/visualizer.json` (вид сверху,
орто-проекция, `parallel_scale` — полу-высота обзора в метрах). Кольцо ~353 м, центр Y≈177.
Для зума на инжекцию — верните `focus/position` к `[7,0.5,*]`, `parallel_scale ~2`.

## 9. Передаваемый интерактивный HTML

`viz_html.py` экспортирует 3D-сцену эксперимента в **самодостаточный `.html`** (plotly встроен в файл):
крутится/зумится в любом браузере **офлайн, без установок**, передаётся как файл. Всё пишется в
**`workspace/results/html/`** (создаётся автоматически). Слои: магниты (цвет по типу), пучок (линии по
частицам), опорный путь (серая ломаная через центры магнитов), BPM (чёрные ромбы на `kind=monitor`).

**Две степени детализации:**
- **детальная** (`_detailed`) — настоящий меш STL, показывает реальную форму магнитов;
- **лёгкая** (`_light`, `--glyph`) — магнит = блок по габаритам STL (12 треуг.), для полного кольца в разы легче.

`--both` пишет **обе** версии сразу (рекомендуется):

```powershell
# ОБЕ версии — лёгкая (блоки) + детальная (STL)
py -3.14 viz_html.py 1   --both --out workspace/results/html/sila_ring.html
#   → sila_ring_light.html (~7 МБ, блоки) + sila_ring_detailed.html (~14 МБ, настоящая форма)
py -3.14 viz_html.py 200 --both --out workspace/results/html/s2lattice.html

# только одна версия
py -3.14 viz_html.py 1   --glyph                     # лёгкая (блоки), быстрый предпросмотр
py -3.14 viz_html.py 200 --decimate 0.5              # детальная с ручным прореживанием
```

| Ключ | Назначение |
|---|---|
| `--both` | Сгенерировать **обе** версии: `<имя>_light.html` (блоки) и `<имя>_detailed.html` (STL). |
| `--glyph` | Только лёгкая: магниты блоками. Для полного кольца (2082 магнита) без прореживания единственный лёгкий вариант. |
| `--decimate 0..1` | Прореживание STL детальной версии (`0.9` = убрать 90%). **`-1` (по умолч.) — авто**: полное кольцо режется под ~300k треуг. (СИЛА: полный STL = 2.3 млн треуг. ⇒ ~190 МБ, авто-режет до 0.87 ≈ 14 МБ, форма сохраняется). **`0` форсит полный STL** (для СИЛА — сотни МБ, только архивом). |
| `--out путь.html` | Своё имя (по умолчанию `html/exp<N>_3d.html`; с `--both` суффиксы `_light`/`_detailed`). |

Заметки: траектория авто-прореживается до ≤8000 точек/частицу; запускать с **`PYTHONIOENCODING=utf-8`**
(иначе финальный лог падает в cp1251-консоли, но HTML уже записан). Детальную версию полного кольца
удобно слать **архивом** (zip ужимает ещё). Пара графиков сравнения ядер `element_bench` пишется
отдельно — `workspace/research/bench/bench_<exp>.html`.

## 10. Замечания

- **`sp_elements.csv` — точка правки токов руками** (шаблон супериода, общие токи → всё
  кольцо). После правки нужен `assemble_lattice_configs.py` (или он вызывается сам в
  `live_visualizer`). ВАЖНО: `assemble` пишет `lattice_configs` без колонок `mat_*` —
  после refresh для матричного ядра прогоните `phys_to_current.py` заново.
- **Согласованность режимов 1/3:** токи должны быть получены из `mat_*` через
  `phys_to_current.py` — тогда диффур (реж.1) и матрица-на-токах (реж.3) сравнимы.
- **Легаси** (скрипты опытов, отчёты, внешние коды) — в `legacy/`. В сценарии не нужны.
- Сборка C++ (при изменении `src/`): CMake → `build/SynchrotronTracker.exe`,
  `SynchrotronOptimizer.exe`, `MatrixTracker.exe`.
