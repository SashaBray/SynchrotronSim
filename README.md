# ⚛️ Synchrotron Particle Tracker v1.0

Relativistic particle tracking engine based on Ansys field maps and 4th order Runge-Kutta integration.
Программный комплекс для релятивистского моделирования траекторий частиц на основе карт полей Ansys и интеграции Рунге-Кутты 4-го порядка.

---

## 📖 Содержание / Contents
1. [Обзор / Overview](#overview)
2. [Файловая структура / File Structure](#file-structure)
3. [Математическая модель / Mathematical Model](#mathematical-model)
4. [Описание таблиц / Table Definitions](#table-definitions)
5. [Сборка и запуск / Build and Run](#build-and-run)

---

## 🔍 Обзор / Overview <a name="overview"></a>

**ENG:** This software is designed for numerical simulation of charged particle trajectories in accelerators. The core feature is the use of real magnetic and electric field maps calculated in Ansys, providing higher accuracy than idealized analytical models.

**RUS:** Программа предназначена для численного моделирования траекторий заряженных частиц в ускорителях. Ключевая особенность — использование реальных карт магнитных и электрических полей, рассчитанных в Ansys, что обеспечивает более высокую точность по сравнению с идеализированными аналитическими моделями.

---

## 📂 Файловая структура / File Structure <a name="file-structure"></a>

*   **src/** — C++ Source code / Исходный код C++.
*   **workspace/** — Working environment / Рабочее пространство.
    *   **tables/** — Configuration CSV tables / Конфигурационные таблицы CSV.
    *   **fields_raw/** — Raw Ansys field maps (.txt) / Исходные карты полей Ansys (.txt).
    *   **fields_bin/** — Processed binary maps for fast loading / Бинарные карты для быстрой загрузки.
    *   **geometry/** — 3D models of devices (.stl) / 3D модели устройств (.stl).
    *   **results/** — Simulation outputs (binary and CSV) / Результаты расчета (бинарные и CSV).

---

## 🧮 Математическая модель / Mathematical Model <a name="mathematical-model"></a>

**ENG:**
*   **Integrator:** 4th order Runge-Kutta method (RK4).
*   **Dynamics:** Relativistic. The particle state is described by the vector $(\vec{x}, \vec{p})$, where $\vec{p}$ is the relativistic momentum. This prevents the particle from exceeding the speed of light.
*   **Forces:** Lorentz force $\vec{F} = q(\vec{E} + \vec{v} \times \vec{B})$ and optional synchrotron radiation energy losses.
*   **Interpolation:** Trilinear interpolation using 8 nodes of a regular grid.
*   **Coordinate System:** $X$ - forward (tangent to orbit), $Y$ - left, $Z$ - up.

**RUS:**
*   **Интегратор:** Метод Рунге-Кутты 4-го порядка (RK4).
*   **Динамика:** Релятивистская. Состояние частицы описывается вектором $(\vec{x}, \vec{p})$, где $\vec{p}$ — релятивистский импульс. Это исключает превышение скорости света.
*   **Силы:** Сила Лоренца $\vec{F} = q(\vec{E} + \vec{v} \times \vec{B})$ и опциональные потери на синхротронное излучение.
*   **Интерполяция:** Трилинейная по 8 узлам регулярной сетки.
*   **Система координат:** $X$ — вперед (касательная к орбите), $Y$ — влево, $Z$ — вверх.

---

## 📊 Описание таблиц / Table Definitions <a name="table-definitions"></a>

**Units of Measurement / Единицы измерения:** SI (Meters, Seconds, Kg, Coulombs, Tesla, Volts) / СИ (Метры, Секунды, Кг, Кулоны, Тесла, Вольты).

### 1. hardware_library.csv

| Column | ENG Description | RUS Описание |
| :--- | :--- | :--- |
| `device_type_id` | Unique name of the device type | Уникальное имя типа устройства |
| `field_type` | Field type (B - magnetic, E - electric) | Тип поля (B - магнитное, E - электр.) |
| `stl_path` | Name of the STL geometry file | Имя файла 3D модели (STL) |
| `nominal_arg` | Current (A) or Voltage (V) of the map | Ток (А) или Напряжение (В) карты |

### 2. lattice_configs.csv

| Column | ENG Description | RUS Описание |
| :--- | :--- | :--- |
| `lattice_id` | ID of the accelerator configuration | ID конфигурации ускорителя |
| `x, y, z` | Device center position [m] | Координаты центра устройства [м] |
| `yaw, pitch, roll` | Rotation angles [rad] | Углы ориентации [рад] |
| `arg_val` | Operational Current (A) or Voltage (V) | Реальный ток (А) или напряжение (В) |

### 3. particle_groups.csv

| Column | ENG Description | RUS Описание |
| :--- | :--- | :--- |
| `group_id` | ID of the particle beam | ID группы (пучка) частиц |
| `m_kg` | Rest mass [kg] | Масса покоя [кг] |
| `q_c` | Charge [C] | Заряд [Кл] |
| `vx, vy, vz` | Initial velocity components [m/s] | Начальная скорость [м/с] |

### 4. experiments.csv

| Column | ENG Description | RUS Описание |
| :--- | :--- | :--- |
| `exp_id` | Unique Experiment ID | Уникальный ID эксперимента |
| `dt` | Integration time step [s] | Шаг интегрирования [с] |
| `max_turns` | Number of turns to simulate | Количество оборотов для расчета |
| `p1..p3 (x,y,z)` | Finish plane definition (3 points) | Координаты финишной плоскости (3 точки) |

---

## 📂 Структура файлов / File Structure

```text
SynchrotronSim/
├── src/                        # C++ Source code (.cpp, .h)
├── workspace/                  # Working directory
│   ├── tables/                 # CSV Configurations
│   ├── fields_raw/             # Raw Ansys .txt maps
│   ├── fields_bin/             # Fast-load binary maps
│   ├── geometry/               # STL models (ignored by git)
│   └── results/                # Simulation outputs
├── visualize.py                # Python 3D Visualizer
└── README.md                   # Documentation
```

---


## 🛠 Сборка и запуск / Build and Run <a name="build-and-run"></a>

**ENG:**
1.  **Build:** Use CMake with MinGW or MSVC.
    ```bash
    mkdir build && cd build
    cmake -G "MinGW Makefiles" ..
    cmake --build .
    ```
2.  **Run:** Execute with experiment ID as argument.
    ```bash
    SynchrotronTracker.exe 1
    ```
3.  **Visualize:** Run the Python script with experiment ID.
    ```bash
    python visualize.py 1
    ```

**RUS:**
1.  **Сборка:** Используйте CMake с MinGW или MSVC.
    ```bash
    mkdir build && cd build
    cmake -G "MinGW Makefiles" ..
    cmake --build .
    ```
2.  **Запуск:** Запустите .exe файл, передав ID эксперимента.
    ```bash
    SynchrotronTracker.exe 1
    ```
3.  **Визуализация:** Запустите Python-скрипт с указанием ID эксперимента.
    ```bash
    python visualize.py 1
    ```

---

## ⚠️ Git Note / Примечание по Git
**ENG:** Binary field maps (`.bin`) and 3D models (`.stl`) are excluded from this repository via `.gitignore` to keep it lightweight. Please provide your own data in the `workspace/` directory.

**RUS:** Бинарные карты полей (`.bin`) и 3D модели (`.stl`) исключены из репозитория через `.gitignore`. Пожалуйста, добавьте собственные данные в папку `workspace/` перед использованием.
