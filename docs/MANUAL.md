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

| Команда | Что делает |
|---|---|
| `py -3.14 lte_to_lattice.py` | Парсит `workspace/lte/*.lte` → `sp_elements.csv` (шаблон супериода) + `lattice_layout.csv` (развёртка кольца). Существующие `arg_val` сохраняются. |
| `py -3.14 assemble_lattice_configs.py` | Разворачивает шаблон по кольцу → `lattice_configs.csv`. **Запускать после любой правки `sp_elements.csv`.** |
| `py -3.14 lte_warmstart.py` | Warm-start физкоэф. `mat_*` из проектных `K` (.lte). |
| `py -3.14 warmstart_ring.py` | Warm-start всего кольца (standard + junction). |
| `py -3.14 build_matrix_db.py` | Сборка вспомогательных таблиц матричной модели. |
| `py -3.14 lte_to_design_orbit.py` | .lte → сегменты опорной орбиты (для оптимизатора). |

## 2. Карты полей

| Команда | Что делает |
|---|---|
| `py -3.14 scan_field_maps.py` | Инвентаризация карт `workspace/fields_bin/`, `fields_raw/`. |
| `py -3.14 gen_synthetic_maps.py` | Генерирует синтетические карты для элементов библиотеки без ANSYS-карты (ресайз эталона под `eff_length`). |

Карты: `.bin` (заголовок 96 б + сетка Bx,By,Bz), `.fld` (текст ANSYS). Скейл поля в трекере:
`map × currentScale × correct_coeff × field_direction`, `currentScale = arg_val / nominal_arg`.

## 3. Сглаживание и симметризация карт

| Команда | Что делает |
|---|---|
| `py -3.14 clean_dipole_maps.py [dev ...]` | Наложить медианную симметрию (Bx,By нечёт / Bz чёт по z), убрать паразитное медианное поле. Без аргументов — все `DP_*` в `synthetic/`. Бэкап `.orig`. |
| `py -3.14 symmetrize_map.py` | Денойз карты по чётностям (симметризация секступоля). *(нужен scipy)* |

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
| `py -3.14 sync_sp_elements.py 101` | Перенести настроенные токи из `lattice_configs.csv` в `sp_elements.csv` (нужно перед `live_visualizer`, т.к. тюн пишет в `lattice_configs`). |
| `py -3.14 reset_magnets.py` | Сброс к дефолту: все `arg_val=0`, `correct_coeff=1`. |
| `py -3.14 manual_kick.py` | Однократный ручной кик токов/коэф. в `sp_elements.csv`. |

**Обзор всего кольца:** камера настроена в `workspace/configs/visualizer.json` (вид сверху,
орто-проекция, `parallel_scale` — полу-высота обзора в метрах). Кольцо ~353 м, центр Y≈177.
Для зума на инжекцию — верните `focus/position` к `[7,0.5,*]`, `parallel_scale ~2`.

## 9. Замечания

- **`sp_elements.csv` — точка правки токов руками** (шаблон супериода, общие токи → всё
  кольцо). После правки нужен `assemble_lattice_configs.py` (или он вызывается сам в
  `live_visualizer`). ВАЖНО: `assemble` пишет `lattice_configs` без колонок `mat_*` —
  после refresh для матричного ядра прогоните `phys_to_current.py` заново.
- **Согласованность режимов 1/3:** токи должны быть получены из `mat_*` через
  `phys_to_current.py` — тогда диффур (реж.1) и матрица-на-токах (реж.3) сравнимы.
- **Легаси** (скрипты опытов, отчёты, внешние коды) — в `legacy/`. В сценарии не нужны.
- Сборка C++ (при изменении `src/`): CMake → `build/SynchrotronTracker.exe`,
  `SynchrotronOptimizer.exe`, `MatrixTracker.exe`.
