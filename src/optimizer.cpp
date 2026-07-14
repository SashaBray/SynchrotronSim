// SynchrotronOptimizer.cpp -- стохастический hill-climber для подбора
// токов магнитов и калибровочных коэффициентов диполей.
//
// Архитектура hot-цикла:
//   - все CSV-файлы читаются ОДИН раз при старте, никакой записи на диск в цикле
//   - field-карты грузятся ОДИН раз (mapCache переиспользуется)
//   - применение параметра -> мутация MagneticElement::currentScale / correctCoeff
//     для заранее посчитанного списка индексов; никакого "assemble"
//   - tracking + loss считаются прямо в памяти (траектория не пишется в bin)
//   - чекпойнты истории в JSON каждые CHECKPOINT_EVERY итераций (на случай Ctrl+C)
//   - финальная JSON-история для viz_optim_history.py
//   - финальные best-значения параметров пишутся обратно в sp_elements.csv и
//     hardware_library.csv (в КОНЦЕ, чтобы visualize.py видел итог)
//   - финальный прогон с best-state перезаписывает trajectories_exp{id}.bin

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <set>
#include <memory>
#include <random>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <filesystem>
#include <system_error>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <tuple>

#ifdef _WIN32
#include <windows.h>
#endif

#include "FieldMap.h"
#include "LatticeManager.h"
#include "Integrator.h"
#include "FinishPlane.h"

namespace fs = std::filesystem;

struct StepTimer {
    std::chrono::high_resolution_clock::time_point t0;
    StepTimer() : t0(std::chrono::high_resolution_clock::now()) {}
    long long ms() const {
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    }
    double sec() const { return ms() / 1000.0; }
};


// =========================================================================
// Загрузка готовой design orbit из CSV (x;y;z), сгенерированного
// prepare_optim_config.py через design_orbit.py. Это ИСТИННАЯ проектная
// траектория (собрана из .lte), а не зигзаг через центры железа магнитов.
// Возвращает пустой вектор при ошибке (вызывающий код может откатиться на
// центры-элементов как fallback).
// =========================================================================
static std::vector<FieldVector> loadDesignOrbitCsv(const std::string& path) {
    std::vector<FieldVector> pts;
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "  design orbit: cannot open " << path << "\n";
        return pts;
    }
    std::string line;
    bool firstRow = true;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (firstRow) {
            firstRow = false;
            // Пропустить заголовок "x;y;z" (или любую первую нечисловую строку)
            if (!line.empty() && !(std::isdigit((unsigned char)line[0]) ||
                                    line[0] == '-' || line[0] == '+' ||
                                    line[0] == '.'))
                continue;
        }
        // Парсим x;y;z (точку с запятой как разделитель)
        double x = 0, y = 0, z = 0;
        size_t p1 = line.find(';');
        if (p1 == std::string::npos) continue;
        size_t p2 = line.find(';', p1 + 1);
        if (p2 == std::string::npos) continue;
        try {
            x = std::stod(line.substr(0, p1));
            y = std::stod(line.substr(p1 + 1, p2 - p1 - 1));
            z = std::stod(line.substr(p2 + 1));
        } catch (...) { continue; }
        pts.push_back({x, y, z});
    }
    return pts;
}

// =========================================================================
// Вход/выход одного магнита (global_id) на эталонной орбите. Используется
// per_element-режимом для построения start_pose и finish_plane на каждом окне.
// =========================================================================
struct MagnetEndpoint {
    int    global_id    = 0;
    std::string device_type_id;
    std::string lte_name;
    double x_in = 0, y_in = 0, yaw_in = 0;
    double x_out = 0, y_out = 0, yaw_out = 0;
};

static std::vector<MagnetEndpoint> loadMagnetEndpointsCsv(const std::string& path) {
    std::vector<MagnetEndpoint> out;
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "  magnet endpoints: cannot open " << path << "\n";
        return out;
    }
    std::string line;
    bool header = true;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (header) {
            // Заголовок начинается с не-цифры
            if (!line.empty() && !std::isdigit((unsigned char)line[0])) {
                header = false;
                continue;
            }
            header = false;
        }
        // Парсим разделители ';'  Колонки:
        // global_id;device_type_id;lte_name;x_in;y_in;yaw_in;x_out;y_out;yaw_out
        std::vector<std::string> cols;
        std::string cur;
        for (char c : line) {
            if (c == ';') { cols.push_back(cur); cur.clear(); }
            else          { cur.push_back(c); }
        }
        cols.push_back(cur);
        if (cols.size() < 9) continue;
        try {
            MagnetEndpoint e;
            e.global_id      = std::stoi(cols[0]);
            e.device_type_id = cols[1];
            e.lte_name       = cols[2];
            e.x_in    = std::stod(cols[3]);
            e.y_in    = std::stod(cols[4]);
            e.yaw_in  = std::stod(cols[5]);
            e.x_out   = std::stod(cols[6]);
            e.y_out   = std::stod(cols[7]);
            e.yaw_out = std::stod(cols[8]);
            out.push_back(std::move(e));
        } catch (...) { continue; }
    }
    return out;
}

// =========================================================================
// Уточнение полилинии аналитическими ДУГАМИ между соседними центрами.
// Порт refine_polyline_with_arcs() из analyze_trajectory_stability.py.
//
// Между p_a (yaw_a) и p_b (yaw_b) при значимом dyaw строится дуга радиуса
// ρ = arc_len / |dyaw|, касательная yaw_a в p_a и yaw_b в p_b. На прямом
// участке (dyaw ≈ 0) — линейная интерполяция. n_subdivide точек на сегмент.
// =========================================================================
static std::vector<FieldVector> refinePolylineWithArcs(
        const std::vector<FieldVector>& positions,
        const std::vector<double>& yaws,
        int n_subdivide)
{
    std::vector<FieldVector> refined;
    if (positions.size() < 2 || n_subdivide < 1 || positions.size() != yaws.size()) {
        return positions;
    }
    refined.reserve(positions.size() * (size_t)n_subdivide + 1);
    refined.push_back(positions[0]);

    const double TWO_PI = 6.283185307179586;
    const double PI     = 3.141592653589793;
    const double DYAW_EPS = 1e-6;
    const double CHORD_EPS = 1e-12;

    for (size_t i = 0; i + 1 < positions.size(); ++i) {
        const FieldVector& p_a = positions[i];
        const FieldVector& p_b = positions[i + 1];
        double yaw_a = yaws[i];
        double yaw_b = yaws[i + 1];
        double dyaw = yaw_b - yaw_a;
        // нормализация dyaw в [-π, π]
        while (dyaw >  PI) dyaw -= TWO_PI;
        while (dyaw < -PI) dyaw += TWO_PI;

        double dx = p_b.x - p_a.x;
        double dy = p_b.y - p_a.y;
        double dz = p_b.z - p_a.z;
        double chord_len = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (std::abs(dyaw) < DYAW_EPS || chord_len < CHORD_EPS) {
            // Прямой сегмент: линейная интерполяция
            for (int k = 1; k <= n_subdivide; ++k) {
                double t = (double)k / n_subdivide;
                refined.push_back({p_a.x + t * dx, p_a.y + t * dy, p_a.z + t * dz});
            }
        } else {
            // Дуга: длина по дуге через формулу для хорды
            double half_dyaw = 0.5 * std::abs(dyaw);
            double arc_len   = chord_len * std::abs(dyaw) / (2.0 * std::sin(half_dyaw));
            double kappa     = dyaw / arc_len;
            double sin_a = std::sin(yaw_a);
            double cos_a = std::cos(yaw_a);
            for (int k = 1; k <= n_subdivide; ++k) {
                double s = (double)k * arc_len / n_subdivide;
                double yaw_s = yaw_a + kappa * s;
                double x_off = (std::sin(yaw_s) - sin_a) / kappa;
                double y_off = (cos_a - std::cos(yaw_s)) / kappa;
                double z_off = dz * ((double)k / n_subdivide);  // z линейно
                refined.push_back({p_a.x + x_off, p_a.y + y_off, p_a.z + z_off});
            }
        }
    }
    return refined;
}

// =========================================================================
// Минимальное расстояние от точки до полилинии (ломаная) -- для оценки
// отклонения частицы от дизайн-орбиты в envelope_loss.
//
// Алгоритм: перебираем сегменты [p_i, p_{i+1}], проекция точки на каждый
// с clamp к [0,1], квадрат расстояния -- берём минимум.
// =========================================================================
static double distanceToPolyline(const FieldVector& p,
                                  const std::vector<FieldVector>& poly)
{
    if (poly.size() < 2) return 0.0;
    double best2 = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i + 1 < poly.size(); ++i) {
        const auto& a = poly[i];
        const auto& b = poly[i+1];
        double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
        double len2 = ux*ux + uy*uy + uz*uz;
        double t = 0.0;
        if (len2 > 1e-30) {
            t = ((p.x - a.x)*ux + (p.y - a.y)*uy + (p.z - a.z)*uz) / len2;
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;
        }
        double cx = a.x + t*ux, cy = a.y + t*uy, cz = a.z + t*uz;
        double dx = p.x - cx, dy = p.y - cy, dz = p.z - cz;
        double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < best2) best2 = d2;
    }
    return std::sqrt(best2);
}

// ============================================================================
// CSV-утилиты (расширены чтением заголовка по сравнению с main.cpp::readCSV).
// ============================================================================
static std::vector<std::string> splitDelim(const std::string& line, char delim) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string v;
    while (std::getline(ss, v, delim)) {
        v.erase(0, v.find_first_not_of(" \t\r\n\""));
        size_t end = v.find_last_not_of(" \t\r\n\"");
        v = (end == std::string::npos) ? "" : v.substr(0, end + 1);
        out.push_back(v);
    }
    return out;
}
static void normalizeEUDecimal(std::string& s) {
    for (char &c : s) if (c == ',') c = '.';
}

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    int col(const std::string& name) const {
        for (int i = 0; i < (int)header.size(); ++i) if (header[i] == name) return i;
        return -1;
    }
};

static CsvTable readCsvTable(const std::string& path) {
    CsvTable t;
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "ERROR: cannot open " << path << std::endl;
        return t;
    }
    std::string line;
    if (!std::getline(f, line)) return t;
    t.header = splitDelim(line, ';');
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto cells = splitDelim(line, ';');
        for (auto& c : cells) normalizeEUDecimal(c);
        t.rows.push_back(std::move(cells));
    }
    return t;
}

// Запись CSV (для финального сохранения best-значений).
static void writeCsvTable(const std::string& path, const CsvTable& t, char delim = ';') {
    std::ofstream f(path);
    if (!f.is_open()) {
        std::cerr << "ERROR: cannot write " << path << std::endl;
        return;
    }
    for (size_t i = 0; i < t.header.size(); ++i) {
        f << t.header[i];
        if (i + 1 < t.header.size()) f << delim;
    }
    f << "\n";
    for (const auto& r : t.rows) {
        for (size_t i = 0; i < r.size(); ++i) {
            f << r[i];
            if (i + 1 < r.size()) f << delim;
        }
        f << "\n";
    }
}

// ============================================================================
// findWorkspace -- скопировано из main.cpp.
// ============================================================================
static std::string findWorkspace(const char* argv0) {
    std::error_code ec;
    std::vector<fs::path> starts;
    #ifdef _WIN32
    char exePath[MAX_PATH];
    DWORD len = GetModuleFileNameA(NULL, exePath, MAX_PATH);
    if (len > 0 && len < MAX_PATH) starts.push_back(fs::path(exePath).parent_path());
    #endif
    if (argv0 && *argv0) {
        fs::path p(argv0);
        if (p.has_parent_path()) starts.push_back(p.parent_path());
    }
    starts.push_back(fs::current_path(ec));
    for (const auto& start : starts) {
        fs::path dir = fs::absolute(start, ec);
        for (int i = 0; i < 10; ++i) {
            fs::path candidate = dir / "workspace";
            if (fs::is_directory(candidate, ec)) return candidate.generic_string() + "/";
            fs::path parent = dir.parent_path();
            if (parent.empty() || parent == dir) break;
            dir = parent;
        }
    }
    return "workspace/";
}

// ============================================================================
// Flat-конфиг (вывод prepare_optim_config.py).
// ============================================================================
struct Matcher {
    std::string column;             // sp_type / sp_element_id / device_type_id
    std::vector<std::string> values; // OR-набор (любое из значений матчит)
};
struct GroupConfig {
    std::string name;
    double stepAbs;
    double stepPercent;
    // target: "CURRENT"     -- меняем arg_val (ток) в sp_elements; el.currentScale
    //         "COEFF_HW"    -- меняем correct_coeff в hardware_library; общий для типа
    //         "COEFF_INST"  -- меняем correct_coeff в sp_elements; per-position
    std::string target;
    std::string file;       // sp_elements.csv / hardware_library.csv
    std::string valueCol;
    std::vector<Matcher> matchers;
};
struct OptimConfig {
    int expId = 0;
    int nIter = 0;
    int seed  = 0;
    int kickEvery = 0;
    double kickScale = 0.0;
    int revertAfterKicks = 3;   // 0 = disable revert
    std::string lossType = "orbit_avg";
    bool lossSquared = true;
    // lossMaxWeight: если > 0, к position-loss добавляется штраф за МАКСИМАЛЬНОЕ
    // отклонение: loss += lossMaxWeight × max_(частица,магнит)(d²) (или d, если
    // lossSquared=false). Заставляет оптимизатор давить худшую точку (пик), а не
    // только сумму. 0 (default) -- выключено.
    double lossMaxWeight = 0.0;
    // interparticleWeight: если > 0, к loss добавляется СУММА КВАДРАТОВ расстояний
    // МЕЖДУ частицами в точках ближайшего подхода к каждой контрольной точке:
    // loss += w · Σ_магнит Σ_{i<j} |p_i - p_j|². Штрафует расхождение пучка
    // (не даёт одной частице отделиться и «утянуть» сумму квадратов pos).
    // Масштаб ~как у pos (тоже по контрольным точкам). 0 (default) -- выключено.
    double interparticleWeight = 0.0;
    // envelopeDecayWeight: если > 0, к loss добавляется штраф за НАКЛОН ТРЕНДА амплитуды
    // бетатронных колебаний. В каждой контрольной точке m: a_m = sqrt(Σ_{i<j}|p_i-p_j|²)
    // -- амплитуда взаимного разлёта частиц. По всем точкам строим ЛИНИЮ ТРЕНДА (МНК-фит
    // a_m ≈ α + β·x, x = m/(M-1) ∈ [0,1]); β -- подъём тренда от начала к концу окна (м).
    // loss += w · max(0, β): штраф ТОЛЬКО за рост (огибающая должна не возрастать, β ≤ 0).
    // Транзиентные пики осцилляций усредняются фитом -- ловится чистый ТРЕНД роста, а не
    // мгновенные превышения. Прим.: симплектичная решётка -- буквального затухания нет,
    // достижим нулевой/отрицательный наклон (стабильная огибающая). 0 (default) -- выкл.
    double envelopeDecayWeight = 0.0;
    // parallelismWeight: если > 0, к loss добавляется штраф за УГЛОВУЮ расходимость
    // пучка: Σ_{i<j}|p̂_i - p̂_j|² (спред единичных направлений импульса частиц на
    // финише окна). Частицы должны двигаться ПАРАЛЛЕЛЬНО. Дополняет позиционный
    // член (тот держит координаты, этот — углы). Активен в ОКОННОМ режиме.
    // 0 (default) — выключено.
    double parallelismWeight = 0.0;
    // svdAngleWeight: вес ПЕРИОДИЧЕСКОГО члена в невязке SVD — направление импульса
    // каждой частицы на финише должно быть ∥ касательной орбиты (fpDir). Только для
    // режима svd (добавляет 2·N_particles записей в residual). 0 = выкл (чистая позиция).
    double svdAngleWeight = 0.0;
    // globalSearchMode: МЕТОД поиска внутри ГЛОБАЛЬНОГО (global) режима.
    //   "monte_carlo" -- случайные greedy-пертурбации групп (метод Монте-Карло;
    //                    алиас "random").
    //   "line_search" -- покоординатный спуск с детерминированным coarse-to-fine
    //                    1D-поиском по диапазону каждого элемента, по ВСЕЙ конфигурации
    //                    (тот же метод последовательного приближения, что в окне, но
    //                    loss глобальный). Параметры поиска (K, max_rounds, диапазоны
    //                    тока/коэффициента) переиспользуются из perElement.* (дефолты
    //                    3/12/100/3, если per_element-блока нет). dynamic_loss и кики
    //                    в этом режиме не применяются.
    // ДЕФОЛТ: "line_search" -- основной метод (детерминированное последовательное
    // приближение по плотному опорному пути). "random" (случайные пертурбации) --
    // опциональная добавка (нужно явно задать global_search_mode="random").
    std::string globalSearchMode = "line_search";
    // Источник координат для position-loss магнитов:
    //   "iron"  (default) -- центры железа из lattice_layout.csv (старое поведение).
    //                        Для не-bending магнитов (KQUAD/KSEXT/KOCT) это же точки
    //                        на эталонной орбите. Для диполей -- центры железа, которые
    //                        в дугах смещены от орбиты на sagitta.
    //   "orbit" -- середина (entry+exit)/2 из magnet_endpoints CSV (точки на
    //              эталонной орбите). Важно если диполи входят в loss.
    std::string positionReference = "iron";
    std::vector<int> lossSpNumbers;
    // Если > 0, в lattice включаем только элементы lattice_layout с
    // global_id <= latticeMaxGid (для экспериментов с "обрезанной" lattice,
    // например пролёт через первые 3 SP с последующим свободным дрейфом).
    int latticeMaxGid = 0;
    // Если > 0, моделируем только первые maxParticles частиц группы (первая —
    // центральная P1 с нулевым смещением). Полезно для фазы настройки поворотных:
    // им нужна только центральная орбита, спутники (для фокусировки линзами)
    // лишь добавляют шум и замедляют. 0 (default) — все частицы группы.
    int maxParticles = 0;
    // device_type_id, исключаемые из набора точек-магнитов для loss
    // (например, постоянные магниты DP_L_PM/DP_sh_PM, которые оптимизатор
    // не подстраивает -- незачем штрафовать орбиту по их центрам).
    std::vector<std::string> lossExcludeDevices;
    // Дополнительные точки-цели для loss (например, центр первого магнита
    // следующего супериода). Задаются строками "TARGET <x> <y> <z> [weight]"
    // в optimconfig. weight по умолчанию 1.0; больше = сильнее штраф этой
    // точки в loss относительно обычных магнитов.
    std::vector<FieldVector> extraTargets;
    std::vector<double>      extraTargetWeights;

    // Точки ПЛОТНОГО опорного пути для composite_next (loss_mode). Задаются
    // строками "PATHPT <gid> <x> <y> <z>" -- центры элементов + уплотнённые
    // прямые участки (вкл. 50м-дрейф). В composite_next loss считает откл. от
    // ЭТИХ точек (а не от центров магнитов), оконно по gid. gid тэгирует, в
    // какое окно точка входит (gid <= gid_last).
    std::vector<FieldVector> pathPoints;
    std::vector<int>         pathPointGids;
    std::vector<double>      pathPointWeights;

    // ============== СОСТАВНОЙ LOSS: угловая компонента ==============
    // Дополнительный штраф за расхождение направления импульса частицы в
    // конце трассы с эталонным направлением (= direction опорной орбиты).
    // Loss_dir = dirWeight × Σ_по_частицам sin²(θ_i),
    // где θ_i -- угол между нормированным импульсом частицы i после
    // финиш-плоскости и dirReference. sin²(θ)=1-(p̂·ref)² для малых θ ≈ θ².
    //
    // Использование: задать SCALAR DIRECTION_REF dx dy dz и
    //                       SCALAR DIRECTION_WEIGHT w
    // в optimconfig. Если dirWeight==0, компонента отключена.
    FieldVector dirReference   = {1.0, 0.0, 0.0};
    double      dirWeight      = 0.0;  // 0 = выключено

    // ============== СОСТАВНОЙ LOSS: огибающая колебаний ==============
    // Штраф за РОСТ амплитуды бетатронных колебаний.
    // Для каждой спутниковой частицы i (P2..P5) сравниваем расстояние до
    // референсной частицы в начале и в конце трекинга:
    //   A_init_i  = |beam0[i].pos - beam0[0].pos|       (~1 мм по построению)
    //   A_final_i = |beam[i].pos  - beam[0].pos|        (после трекинга)
    //   growth_i  = max(0, A_final_i - A_init_i)
    // Loss_env = envWeight × Σ_i growth_i².
    // Затухание (A_final < A_init) НЕ наказывается -- это полезный эффект.
    double envWeight = 0.0;

    // ============== СПОСОБ ПОСТРОЕНИЯ design-полилинии ==============
    // "linear" (default): прямая ломаная через центры всех элементов.
    //     Между двумя центрами с разным yaw оказывается ХОРДА, а не дуга
    //     дипольного магнита — sagitta-артефакт.
    // "arcs": между парами соседних центров строится АНАЛИТИЧЕСКАЯ ДУГА
    //     по yaw_i, yaw_{i+1} (kappa = (yaw_{i+1} - yaw_i) / arc_len).
    //     Это физически корректная опорная орбита; envelope_loss
    //     становится «честным» — расстояние до проектной кривой, а не
    //     до её хордовой аппроксимации.
    // envArcSubdivide — число подразбиений на каждый сегмент при "arcs".
    std::string envPolylineMode = "linear";
    int envArcSubdivide = 10;

    // ============== Истинная design orbit из .lte ==============
    // Путь к CSV (x;y;z) с заранее собранной проектной траекторией. Генерируется
    // Python-скриптом prepare_optim_config.py через модуль design_orbit.py из БД
    // сегментов (workspace/tables/design_orbit_segments.csv) -- которая, в свою
    // очередь, получена парсингом .lte (lte_to_design_orbit.py). Это физически
    // корректная орбита (не зависит от геометрии железа в lattice_configs).
    //
    // Если задано и файл читается -- используется ВМЕСТО построения через
    // центры элементов (envPolylineMode игнорируется).
    std::string designOrbitCsv = "";

    std::string historyOut = "optim_history.json";
    std::vector<GroupConfig> groups;

    // ============== РЕЖИМ LOSS-ФУНКЦИИ ==============
    // "composite" (default): position + direction + envelope.
    //     Полная составная loss для финальной полировки.
    // "orbit": только position (sum d² по магнитам опорной орбиты).
    //     Используется для грубой первичной настройки, когда мы хотим
    //     максимально быстро привести пучок к опорной орбите без оглядки
    //     на угол выхода и рост бетатронных колебаний.
    // "envelope": только envelope (рост амплитуды бетатронных колебаний).
    //     Используется когда орбита уже хороша, но колебания спутников
    //     надо подавить — оптимизатор будет настраивать линзы под условие
    //     A_final ≤ A_init по всем спутникам.
    // CLI-флаг --loss=orbit|composite|envelope перебивает JSON.
    // Реализуется зануление dirWeight/envWeight для "orbit" и зануление
    // usePositionLoss/dirWeight для "envelope".
    std::string lossMode = "composite";
    // Если false — position-loss (sum d² по магнитам и TARGET'ам) не
    // включается в outLoss. Активируется при lossMode == "envelope".
    bool usePositionLoss = true;

    // ============== РЕЖИМ ОПТИМИЗАЦИИ ==============
    // method: РЕЖИМ настройки.
    //   "global"   (default) -- глобальный: настраивается вся конфигурация сразу.
    //                            Метод внутри режима задаёт globalSearchMode
    //                            (line_search | monte_carlo).
    //   "windowed"            -- оконный: скользящее окно по магнитам. Метод внутри
    //                            окна задаёт perElement.searchMode.
    //   "gradient"            -- градиентный (FD + Adam).
    // Переопределяется CLI-флагом --method=/--mode=. Историч. алиасы montecarlo/mc
    // -> global, per_element/pe -> windowed принимаются (см. canonicalMode).
    std::string method = "global";

    // ============== DYNAMIC LOSS REWEIGHTING (только MC) ==============
    // Веса компонентов composite (position, direction, envelope) вращаются по
    // синусоиде с равномерным распределением фаз. Идея: в каждый момент
    // ВРЕМЕНИ один компонент имеет повышенный вес, другой -- пониженный.
    // Это позволяет MC выбираться из локальных минимумов, в которых застревает
    // статически-взвешенная composite loss.
    //
    // Формула:
    //   phase_i = 2π × i / 3                        (i = 0=pos, 1=dir, 2=env)
    //   w_i(t) = base_i × [1 + A × cos(2π·t/T - phase_i)]
    //   L_dyn(t) = w_pos(t)·pos + w_dir(t)·dir + w_env(t)·env
    //
    // MC принимает шаг по L_dyn (новый < текущего с теми же весами времени t).
    // В history.bestLoss пишется СТАТИЧЕСКИЙ composite (A=0) -- для объективного
    // отслеживания прогресса.
    //
    // periodIterations: если 0 -> авто 5×N_params (каждый параметр в среднем
    // мутирует ~5 раз за полный оборот фокуса).
    struct DynamicLossConfig {
        bool   enabled         = false;
        double amplitude       = 0.7;   // A ∈ [0, 1]
        int    periodIterations = 0;    // T (0 = auto)
    } dynamicLoss;

    // ============== ОКОННЫЙ РЕЖИМ (windowed, скользящее окно) ==============
    // Активен при method == "windowed". Окно из W магнитов скользит по
    // global_id (от первого активного до последнего). На каждом положении
    // окна — m итераций hill-climber MC, варьируя только группы, чьи
    // elementIndex.global_id ∈ окно. Между окнами найденные значения
    // сохраняются (persistent state). N_passes таких проходов по решётке.
    struct PerElementConfig {
        int windowSize     = 5;   // W: магнитов в окне
        int windowStride   = 1;   // сдвиг окна за один шаг (1 = максимальное перекрытие)
        int itersPerWindow = 200; // m: MC-итераций в каждом положении окна
        int passes         = 1;   // N: сколько раз пройти всю решётку
        // Если true -- частицы стартуют с начала СУПЕРИОДА, в котором находится
        // первый магнит окна (мы прокачиваем пучок через все предыдущие магниты SP).
        // Учитывает накопленные ошибки и физически соответствует "реальному" пучку.
        // Если false (default) -- идеальный старт сразу перед окном на эталонной орбите.
        bool startFromSpStart = false;

        // Штраф за отклонение P1 ПО ДОРОГЕ от design_polyline (а не только в
        // финиш-плоскости окна). Решает проблему "пучок срезает дугу диполей":
        // без него оптимизатор может прийти к решению где P1 летит почти прямо
        // через дугу и попадает в target с малой ошибкой, физически отклоняясь
        // на 5+ см в середине дуги.
        //   path_loss = Σ_steps d²(P1, design_polyline)  (sample каждый
        //               path_sample_stride-й шаг пушера)
        // Режимы:
        //   pathOnly=false (default): legacy. Loss = ‖P1−T‖² + 0.25·sat + dirW·dir + pathW·path
        //   pathOnly=true: Loss = path_loss (path_weight игнорируется -- считается = 1).
        //                  «Чистая сумма квадратов отклонений от опорной орбиты».
        bool   pathOnly         = false;
        double pathWeight       = 0.0;
        int    pathSampleStride = 10;

        // Режим loss для per-window simulation:
        //   "endpoint" (default)  -- ‖P1−T‖² + 0.25·sat + dirW·dir + pathW·path
        //   "path_only"           -- Σ d²(P1, design_polyline)
        //   "composite"           -- composite loss для магнитов В ОКНЕ.
        //                             Симуляция через окно (с mask).
        //   "composite_full"      -- composite loss с учётом ВСЕГО пути от 0
        //                             до финиш-плоскости окна. mask отключена,
        //                             magnetPos = все lattice-магниты с gid ≤
        //                             gid_last окна. Учитывает ошибки прохода
        //                             пучка ДО окна.
        std::string lossMode = "endpoint";

        // accumulateGroups: если true, набор активных групп РАСТЁТ внутри
        // одного pass'а. Окно 1 видит свои ручки, окно 2 -- свои + ручки окна 1,
        // последнее окно видит все. Между pass'ами набор сбрасывается.
        // По умолчанию false (стандартный sliding window).
        bool accumulateGroups = false;

        // itersPerParam: если > 0, число итераций hill-climber для окна =
        // itersPerParam × N_active_groups. Полезно при accumulateGroups=true,
        // когда число параметров растёт. Если 0 (default) -- используется
        // фиксированное itersPerWindow.
        int  itersPerParam   = 0;

        // totalIters: если > 0, ОБЩИЙ бюджет симуляций на проход распределяется
        // РАВНОМЕРНО по окнам: m_iter каждого окна = totalIters / num_windows.
        // Это убирает геометрический рост к концу супериода (когда при
        // accumulate N_active растёт и itersPerParam×N_active взрывается).
        // Имеет приоритет над itersPerParam. 0 (default) -- выключено.
        int  totalIters      = 0;

        // itersPerNew: если > 0 (вместе с totalIters), распределение бюджета по
        // окнам "новый магнит — фиксированно, старые — по остатку":
        //   m_iter(окно) = N_new × itersPerNew + остаток/num_windows,
        //   остаток = max(0, totalIters − G×itersPerNew), G = всего групп.
        // Новые элементы окна получают гарантированный (большой) бюджет, уже
        // настроенные — делят остаток. Приоритет над простым равномерным режимом.
        int  itersPerNew     = 0;

        // groupDipoleArrays: если true, при подходе окна к ЛИНЕЙКЕ диполей
        // (подряд идущие dipole-slice'ы, device_type начинается с "DP") окно
        // сразу захватывает ВСЮ линейку целиком -- слайсы одной системы входят
        // в настройку совместно (одним окном), а не по одному. Внутренние окна
        // линейки пропускаются (не тратим бюджет на по-слайсовую настройку).
        bool groupDipoleArrays = false;

        // advanceMeanThresholdMm: гейт продвижения окна. Если > 0, окно НЕ
        // сдвигается дальше, пока среднее отклонение траекторий частиц от
        // опорной траектории (mean, мм) не упадёт <= порога. Сначала
        // выполняется обычный бюджет (m_iter_dynamic), затем — если mean всё
        // ещё выше порога — окно продолжает оптимизироваться до достижения
        // порога ИЛИ до защитного потолка (advanceMaxIterMult × m_iter_dynamic).
        // Решает проблему: малый шаг диполей не успевает проложить траекторию
        // за фиксированный бюджет. 0 (default) -- гейт выключен.
        // Mean для гейта считается по lattice-магнитам в области симуляции
        // (gid > 0), БЕЗ extra TARGET'ов на дрейфе.
        double advanceMeanThresholdMm = 0.0;
        int    advanceMaxIterMult     = 10;  // защитный потолок = mult × m_iter

        // lensLagWindows: режим "lagged-lens". Если > 0, группы с device_type из
        // lensLaggedDeviceTypes (линзы/фокусирующие элементы) включаются в набор
        // тюнинга НЕ сразу при входе их магнита в окно, а на lensLagWindows окон
        // позже — когда поворотные элементы в этой области уже устаканились.
        // Реализует двухфазность "поворотные → линзы" локально, в один проход:
        // фронт окна кладёт траекторию поворотными, а в хвосте (на lag окон
        // позади) включаются линзы на фокусировку уже-проложенной орбиты.
        // 0 (default) — выключено (все группы активируются сразу, как обычно).
        int lensLagWindows = 0;
        // device_type_id, считающиеся "линзами" (отложенная активация).
        std::set<std::string> lensLaggedDeviceTypes;

        // ---- Режим поиска внутри окна ----
        // searchMode: "random" (default) — случайные greedy-пертурбации групп
        //             (исторический режим). "line_search" — покоординатный спуск
        //             с детерминированным coarse-to-fine 1D-поиском по каждому
        //             элементу (метод последовательного приближения): K точек
        //             равномерно по диапазону, зум вокруг лучшей точки (включая
        //             incumbent), уполовинивание интервала; n раундов на элемент
        //             распределяется из бюджета окна равномерно, остаток — на
        //             повторные свипы по группам.
        // ДЕФОЛТ: "line_search" -- основной метод. "random" (свипы/пертурбации) --
        // опциональная добавка (нужно явно задать search_mode="random"). Кики тоже
        // опциональны (kick_every>0).
        std::string searchMode = "line_search";
        int    lineSearchPoints   = 3;     // K точек на раунд (>=3 для сходимости)
        int    lineSearchMaxRounds = 12;   // потолок раундов на элемент за свип
        int    lineSearchMinRounds = 1;    // ГАРАНТИРОВАННЫЙ нижний порог раундов на элемент за визит
        double currentRangeAmp    = 100.0; // |ток| <= 100 А (диапазон [-100,100])
        double coeffRange         = 3.0;   // |coeff| <= 3 (диапазон [-3,3])
        // ЛОКАЛЬНЫЙ радиус поиска: если > 0, line-search ищет НЕ во всём диапазоне
        // [-range,range], а в ОКРЕСТНОСТИ вокруг текущего значения: [cur-radius,
        // cur+radius] (зажатой в [-range,range]). Лекарство от railing'а: токи
        // двигаются инкрементально от текущего, а не прыгают на упор. За свипы
        // окрестность "переезжает" (центрируется на обновлённом значении), так что
        // большой ток достижим постепенно. 0 = искать во всём диапазоне (как раньше).
        double lineSearchLocalCurrent = 0.0;  // ± радиус для токов (А)
        double lineSearchLocalCoeff   = 0.0;  // ± радиус для коэффициентов
        // КИК в глобальном line-search: каждые kickEverySweeps свипов встряхнуть ВСЕ
        // группы на ±kickAmplitude×|значение| (доля; 0.01 = 1%), принять безусловно,
        // затем следующий свип доводит из встряхнутого состояния. bestValues хранит
        // лучшее, поэтому плохой кик не теряет результат. 0 = выкл (как раньше).
        int    lineSearchKickEverySweeps = 0;
        double lineSearchKickAmplitude   = 0.0;
        // БЛОЧНЫЙ покоординатный спуск: blockSize групп тюнятся ВМЕСТЕ за раз (по
        // очереди блоков). На блок -- оценка градиента (конечн. разности) + line-search
        // вдоль -grad. Двигает n связанных токов совместно -> выходит из координатных
        // застреваний, недоступных по одному. 1 = по одной группе (как раньше).
        int    lineSearchBlockSize = 1;
        // ОТКАТ (accept-final): по умолчанию финал = BEST (монотонно, откат невозможен).
        // true -> сохраняется ТЕКУЩЕЕ исследованное состояние (даже хуже best), чтобы
        // система могла «шагнуть назад ради двух вперёд». Безопасный кап: если текущее
        // катастрофически разлетелось (loss >> best), всё равно откат к best.
        bool   lineSearchAcceptFinal = false;
        // FIRST-IMPROVEMENT: в global line-search принимать ПЕРВОЕ значение группы, улучшившее
        // тюн, и не дозумивать (экономит бюджет + грубее => реже застревает в остром экстремуме).
        // false (default) -> полный coarse-to-fine зум к локальному минимуму группы.
        bool   lineSearchFirstImprove = false;
        // СПЛАЙН-ШАГ: вместо coarse-to-fine зума — K проб по окрестности, натуральный кубический
        // сплайн loss(value) и его минимум на отрезке как кандидат (сглаживает дискретизац. шум,
        // находит минимум МЕЖДУ пробами). Кандидат всё равно проверяется реальной симуляцией.
        bool   lineSearchSpline = false;
        // ПАРАБОЛА (МНК по всем K пробам): least-squares y≈Au²+Bu+C, вершина = кандидат.
        // Защита: только если парабола ОТКРЫТА ВВЕРХ (минимум) и вершина в диапазоне проб —
        // иначе откат на лучшую пробу. Усредняет шум, не ныряет в минус как кубич. сплайн.
        bool   lineSearchParabola = false;
        // ---- Режим LENS (оконный по линзам, свежий прямой пучок) ----
        // windowUnit: "magnet" (скользящее окно по всем активным магнитам, как
        // раньше) или "lens" (окно по ЛИНЗАМ; в начало окна инжектируется свежий
        // прямой несжатый пучок). lensWindowing включается при method=="lens".
        std::string windowUnit = "magnet";
        // Типы устройств, считающиеся ЛИНЗАМИ (префикс-матч, как isDipGid по "DP").
        // Пусто => дефолт {"QP","DQ"}. Окно содержит lensWindowSize таких линз.
        std::set<std::string> lensDeviceTypes;
        int    lensWindowSize = 4;       // N линз в окне (опционально)
        double injectSpreadM  = 0.001;   // ± поперечный спред инжектируемого пучка (м)
        // ---- Прогрессивное расширение в дрейф ----
        // driftExtensionM: если > 0, ПОСЛЕ оконной укладки SP1 запускается фаза
        // "вползания" в дрейф: рубеж растёт по driftChunkM метров (2, 4, ... до
        // driftExtensionM). На каждом рубеже контрольные точки = весь путь до рубежа
        // (SP1 + дрейф), финиш-плоскость на рубеже; настраиваются ВСЕ группы SP1
        // (line-search, локальный радиус). Пучок обязан следовать линии непрерывно
        // до растущего рубежа -> ни "выстрела", ни расщепления. Точки дрейфа в пути
        // тэгируются прогрессивными gid'ами (см. build_reference_path drift_chunk_m).
        double driftExtensionM    = 0.0;   // до скольки метров тянуть (0 = фаза выкл.)
        double driftChunkM        = 2.0;   // шаг рубежа (м)
        int    driftItersPerChunk = 480;   // бюджет оценок на рубеж (line-search)
    } perElement;

    // Путь к CSV с вход/выход каждого магнита на эталонной орбите. Нужен
    // per_element-режиму для построения start_pose и finish_plane на окне.
    std::string magnetEndpointsCsv = "";

    // Параметры градиентного режима. Активны только при method == "gradient".
    struct GradientConfig {
        std::string fdMethod = "central";   // "central" (точнее, 2M вызовов) | "forward" (M+1)
        int    nIter        = 500;
        double learningRate = 0.05;
        // Правило обновления:
        //   "normalized_gd" (default) — Δθ_j = -lr · g_j / max_k|g_k|.
        //     Шаг по строгому направлению -∇L, нормированному по ∞-норме градиента.
        //     Параметр с максимальным |g| смещается на ±lr, остальные пропорционально.
        //   "adam" — классический Adam с моментами m, v и bias-correction.
        //     На малом бюджете итераций Adam работает хуже: на первой итерации
        //     даёт sign-based шаг ±lr по каждой координате, что для функций с
        //     неравной чувствительностью по θ часто оказывается не вдоль -∇L.
        std::string updateRule = "normalized_gd";
        // Adam-моменты (только для updateRule == "adam")
        double beta1  = 0.9;
        double beta2  = 0.999;
        double epsAdam= 1e-8;
        // ε для конечных разностей; меньше, чем шаг MC (иначе FD выходит из локала).
        double fdStepAbs     = 0.01;   // для arg_val, А
        double fdStepPercent = 0.001;  // для correct_coeff, доля от текущего |θ|
        // Backtracking line search всегда включён. max_halvings контролирует, сколько раз
        // мы делим α пополам, прежде чем считать итерацию неудачной.
        int    maxHalvings = 5;
    } gradient;
};

// ============================================================================
// Канонизация имени РЕЖИМА оптимизации.
//   Режимы: "windowed" (оконный) | "global" (глобальный) | "gradient".
//   Историч. алиасы: montecarlo/mc -> global; per_element/perel/pe -> windowed;
//   grad -> gradient.
// Прим.: "montecarlo" исторически называл ГЛОБАЛЬНЫЙ режим. Теперь Монте-Карло --
// это МЕТОД (случайные кики), задаётся через global_search_mode/search_mode;
// см. canonicalSearchMethod ниже.
// ============================================================================
static std::string canonicalMode(std::string m) {
    for (char& c : m) if (c >= 'A' && c <= 'Z') c += 32;  // lowercase
    if (m == "global" || m == "montecarlo" || m == "mc")              return "global";
    if (m == "windowed" || m == "window" || m == "per_element"
        || m == "perel" || m == "pe")                                return "windowed";
    if (m == "lens" || m == "lenses")                                return "lens";
    if (m == "gradient" || m == "grad")                              return "gradient";
    return m;  // неизвестно -> как есть (вызывающий проверит и ошибётся)
}
// Канонизация имени МЕТОДА поиска внутри режима.
//   "line_search" (последовательное приближение) | "monte_carlo" (случайные кики).
//   Алиасы метода Монте-Карло: random/montecarlo/mc.
static std::string canonicalSearchMethod(std::string m) {
    for (char& c : m) if (c >= 'A' && c <= 'Z') c += 32;  // lowercase
    if (m == "line_search" || m == "linesearch" || m == "ls")        return "line_search";
    if (m == "monte_carlo" || m == "montecarlo" || m == "mc"
        || m == "random")                                            return "monte_carlo";
    if (m == "svd" || m == "response_matrix" || m == "gauss_newton"
        || m == "gn")                                                return "svd";
    return m;
}

static OptimConfig loadFlatConfig(const std::string& path) {
    OptimConfig cfg;
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "ERROR: cannot open " << path << std::endl;
        return cfg;
    }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        // SCALAR <key> <value...>
        if (line.rfind("SCALAR ", 0) == 0) {
            std::string rest = line.substr(7);
            auto sp = rest.find(' ');
            std::string key = (sp == std::string::npos) ? rest : rest.substr(0, sp);
            std::string val = (sp == std::string::npos) ? ""   : rest.substr(sp + 1);
            // Стрипуем whitespace со значения: prepare_optim_config.py выравнивает
            // SCALAR'ы пробелами для читабельности, и без strip'а строковые значения
            // (например "central") приходят как "      central" → сравнение ломается.
            auto lstrip = val.find_first_not_of(" \t");
            if (lstrip == std::string::npos) val = "";
            else                              val = val.substr(lstrip);
            auto rstrip = val.find_last_not_of(" \t\r\n");
            if (rstrip != std::string::npos) val = val.substr(0, rstrip + 1);
            if      (key == "EXP_ID")        cfg.expId      = std::stoi(val);
            else if (key == "N_ITER")        cfg.nIter      = std::stoi(val);
            else if (key == "SEED")          cfg.seed       = std::stoi(val);
            else if (key == "KICK_EVERY")    cfg.kickEvery  = std::stoi(val);
            else if (key == "KICK_SCALE")    cfg.kickScale  = std::stod(val);
            else if (key == "REVERT_AFTER_KICKS") cfg.revertAfterKicks = std::stoi(val);
            else if (key == "LOSS_TYPE")     cfg.lossType   = val;
            else if (key == "LOSS_SQUARED")  cfg.lossSquared= (std::stoi(val) != 0);
            else if (key == "LOSS_MAX_WEIGHT") cfg.lossMaxWeight = std::stod(val);
            else if (key == "INTERPARTICLE_WEIGHT") cfg.interparticleWeight = std::stod(val);
            else if (key == "ENVELOPE_DECAY_WEIGHT") cfg.envelopeDecayWeight = std::stod(val);
            else if (key == "PARALLELISM_WEIGHT")    cfg.parallelismWeight = std::stod(val);
            else if (key == "SVD_ANGLE_WEIGHT")      cfg.svdAngleWeight = std::stod(val);
            else if (key == "GLOBAL_SEARCH_MODE") cfg.globalSearchMode = canonicalSearchMethod(val);
            else if (key == "LOSS_POSITION_REFERENCE") cfg.positionReference = val;
            else if (key == "LOSS_SP") {
                cfg.lossSpNumbers.clear();
                std::stringstream ss(val);
                int v;
                while (ss >> v) cfg.lossSpNumbers.push_back(v);
            }
            else if (key == "LATTICE_MAX_GID") cfg.latticeMaxGid = std::stoi(val);
            else if (key == "MAX_PARTICLES")   cfg.maxParticles  = std::stoi(val);
            else if (key == "DIRECTION_REF") {
                std::stringstream ss(val);
                double dx, dy, dz;
                if (ss >> dx >> dy >> dz) {
                    double n = std::sqrt(dx*dx + dy*dy + dz*dz);
                    if (n > 1e-30) { dx /= n; dy /= n; dz /= n; }
                    cfg.dirReference = {dx, dy, dz};
                }
            }
            else if (key == "DIRECTION_WEIGHT") cfg.dirWeight = std::stod(val);
            else if (key == "ENVELOPE_WEIGHT")  cfg.envWeight = std::stod(val);
            else if (key == "ENV_POLYLINE_MODE") cfg.envPolylineMode = val;
            else if (key == "ENV_ARC_SUBDIVIDE") cfg.envArcSubdivide = std::stoi(val);
            else if (key == "LOSS_EXCLUDE_DEVICES") {
                cfg.lossExcludeDevices.clear();
                std::stringstream ss(val);
                std::string dev;
                while (ss >> dev) cfg.lossExcludeDevices.push_back(dev);
            }
            else if (key == "HISTORY_OUT")   cfg.historyOut = val;
            // -------- loss mode --------
            else if (key == "LOSS_MODE")             cfg.lossMode = val;
            // -------- gradient mode --------
            else if (key == "METHOD")                cfg.method = canonicalMode(val);
            else if (key == "GRAD_FD_METHOD")        cfg.gradient.fdMethod = val;
            else if (key == "GRAD_UPDATE_RULE")      cfg.gradient.updateRule = val;
            else if (key == "GRAD_N_ITER")           cfg.gradient.nIter = std::stoi(val);
            else if (key == "GRAD_LEARNING_RATE")    cfg.gradient.learningRate = std::stod(val);
            else if (key == "GRAD_BETA1")            cfg.gradient.beta1 = std::stod(val);
            else if (key == "GRAD_BETA2")            cfg.gradient.beta2 = std::stod(val);
            else if (key == "GRAD_EPS")              cfg.gradient.epsAdam = std::stod(val);
            else if (key == "GRAD_FD_STEP_ABS")      cfg.gradient.fdStepAbs = std::stod(val);
            else if (key == "GRAD_FD_STEP_PERCENT")  cfg.gradient.fdStepPercent = std::stod(val);
            else if (key == "GRAD_MAX_HALVINGS")     cfg.gradient.maxHalvings = std::stoi(val);
            // -------- dynamic loss reweighting (MC only) --------
            else if (key == "DYN_LOSS_ENABLED")       cfg.dynamicLoss.enabled = (std::stoi(val) != 0);
            else if (key == "DYN_LOSS_AMPLITUDE")     cfg.dynamicLoss.amplitude = std::stod(val);
            else if (key == "DYN_LOSS_PERIOD")        cfg.dynamicLoss.periodIterations = std::stoi(val);
            // -------- per_element mode --------
            else if (key == "PER_ELEM_WINDOW_SIZE")       cfg.perElement.windowSize     = std::stoi(val);
            else if (key == "PER_ELEM_WINDOW_STRIDE")     cfg.perElement.windowStride   = std::stoi(val);
            else if (key == "PER_ELEM_ITERS_PER_WINDOW")  cfg.perElement.itersPerWindow = std::stoi(val);
            else if (key == "PER_ELEM_PASSES")            cfg.perElement.passes         = std::stoi(val);
            else if (key == "PER_ELEM_START_FROM_SP_START") cfg.perElement.startFromSpStart = (std::stoi(val) != 0);
            else if (key == "PER_ELEM_PATH_ONLY")         cfg.perElement.pathOnly = (std::stoi(val) != 0);
            else if (key == "PER_ELEM_PATH_WEIGHT")       cfg.perElement.pathWeight = std::stod(val);
            else if (key == "PER_ELEM_PATH_SAMPLE_STRIDE") cfg.perElement.pathSampleStride = std::stoi(val);
            else if (key == "PER_ELEM_LOSS_MODE")          cfg.perElement.lossMode = val;
            else if (key == "PER_ELEM_ACCUMULATE_GROUPS")  cfg.perElement.accumulateGroups = (std::stoi(val) != 0);
            else if (key == "PER_ELEM_ITERS_PER_PARAM")    cfg.perElement.itersPerParam = std::stoi(val);
            else if (key == "PER_ELEM_TOTAL_ITERS")        cfg.perElement.totalIters = std::stoi(val);
            else if (key == "PER_ELEM_ITERS_PER_NEW")      cfg.perElement.itersPerNew = std::stoi(val);
            else if (key == "PER_ELEM_GROUP_DIPOLE_ARRAYS") cfg.perElement.groupDipoleArrays = (std::stoi(val) != 0);
            else if (key == "PER_ELEM_ADVANCE_MEAN_MM")    cfg.perElement.advanceMeanThresholdMm = std::stod(val);
            else if (key == "PER_ELEM_ADVANCE_MAX_ITER_MULT") cfg.perElement.advanceMaxIterMult = std::stoi(val);
            else if (key == "PER_ELEM_LENS_LAG_WINDOWS")   cfg.perElement.lensLagWindows = std::stoi(val);
            else if (key == "PER_ELEM_SEARCH_MODE")        cfg.perElement.searchMode = canonicalSearchMethod(val);
            else if (key == "PER_ELEM_LS_POINTS")          cfg.perElement.lineSearchPoints = std::stoi(val);
            else if (key == "PER_ELEM_LS_MAX_ROUNDS")      cfg.perElement.lineSearchMaxRounds = std::stoi(val);
            else if (key == "PER_ELEM_LS_MIN_ROUNDS")      cfg.perElement.lineSearchMinRounds = std::stoi(val);
            else if (key == "PER_ELEM_CURRENT_RANGE")      cfg.perElement.currentRangeAmp = std::stod(val);
            else if (key == "PER_ELEM_COEFF_RANGE")        cfg.perElement.coeffRange = std::stod(val);
            else if (key == "PER_ELEM_LS_LOCAL_CURRENT")   cfg.perElement.lineSearchLocalCurrent = std::stod(val);
            else if (key == "PER_ELEM_LS_LOCAL_COEFF")     cfg.perElement.lineSearchLocalCoeff = std::stod(val);
            else if (key == "PER_ELEM_LS_KICK_EVERY_SWEEPS") cfg.perElement.lineSearchKickEverySweeps = std::stoi(val);
            else if (key == "PER_ELEM_LS_KICK_AMPLITUDE")  cfg.perElement.lineSearchKickAmplitude = std::stod(val);
            else if (key == "PER_ELEM_LS_BLOCK_SIZE")      cfg.perElement.lineSearchBlockSize = std::stoi(val);
            else if (key == "PER_ELEM_LS_ACCEPT_FINAL")    cfg.perElement.lineSearchAcceptFinal = (std::stoi(val) != 0);
            else if (key == "PER_ELEM_LS_FIRST_IMPROVE")   cfg.perElement.lineSearchFirstImprove = (std::stoi(val) != 0);
            else if (key == "PER_ELEM_LS_SPLINE")          cfg.perElement.lineSearchSpline = (std::stoi(val) != 0);
            else if (key == "PER_ELEM_LS_PARABOLA")        cfg.perElement.lineSearchParabola = (std::stoi(val) != 0);
            else if (key == "PER_ELEM_DRIFT_EXTENSION_M")  cfg.perElement.driftExtensionM = std::stod(val);
            else if (key == "PER_ELEM_DRIFT_CHUNK_M")      cfg.perElement.driftChunkM = std::stod(val);
            else if (key == "PER_ELEM_DRIFT_ITERS_PER_CHUNK") cfg.perElement.driftItersPerChunk = std::stoi(val);
            else if (key == "PER_ELEM_WINDOW_UNIT")        cfg.perElement.windowUnit = val;
            else if (key == "PER_ELEM_LENS_WINDOW_SIZE")   cfg.perElement.lensWindowSize = std::stoi(val);
            else if (key == "PER_ELEM_INJECT_SPREAD_M")    cfg.perElement.injectSpreadM = std::stod(val);
            continue;
        }
        // STRING <key> <value with possible spaces / slashes>
        // Используется для строковых параметров, которые могут содержать
        // пробелы или Windows-пути (в SCALAR пробелы являются разделителями).
        if (line.rfind("STRING ", 0) == 0) {
            std::string rest = line.substr(7);
            auto sp = rest.find(' ');
            if (sp == std::string::npos) continue;
            std::string key = rest.substr(0, sp);
            std::string val = rest.substr(sp + 1);
            auto lstrip = val.find_first_not_of(" \t");
            if (lstrip == std::string::npos) val = "";
            else                              val = val.substr(lstrip);
            auto rstrip = val.find_last_not_of(" \t\r\n");
            if (rstrip != std::string::npos) val = val.substr(0, rstrip + 1);
            if (key == "DESIGN_ORBIT_CSV")     cfg.designOrbitCsv = val;
            else if (key == "MAGNET_ENDPOINTS_CSV") cfg.magnetEndpointsCsv = val;
            else if (key == "HISTORY_OUT")     cfg.historyOut = val;
            else if (key == "LENS_LAGGED_DEVICE_TYPES") {
                // comma-separated device_type_id, считающиеся линзами (lagged).
                std::stringstream ss(val);
                std::string tok;
                while (std::getline(ss, tok, ',')) {
                    auto a = tok.find_first_not_of(" \t");
                    auto b = tok.find_last_not_of(" \t\r\n");
                    if (a != std::string::npos)
                        cfg.perElement.lensLaggedDeviceTypes.insert(tok.substr(a, b - a + 1));
                }
            }
            else if (key == "LENS_DEVICE_TYPES") {
                // comma-separated префиксы device_type_id, считающиеся ЛИНЗАМИ
                // (границы окон в режиме lens). Префикс-матч (QP -> QP_L_1A и т.п.).
                std::stringstream ss(val);
                std::string tok;
                while (std::getline(ss, tok, ',')) {
                    auto a = tok.find_first_not_of(" \t");
                    auto b = tok.find_last_not_of(" \t\r\n");
                    if (a != std::string::npos)
                        cfg.perElement.lensDeviceTypes.insert(tok.substr(a, b - a + 1));
                }
            }
            continue;
        }
        // TARGET <x> <y> <z> [weight]  -- доп. точка-магнит для loss
        // weight (4-й аргумент) -- множитель вклада этой точки. По умолчанию 1.0.
        if (line.rfind("TARGET ", 0) == 0) {
            std::stringstream ss(line.substr(7));
            double x, y, z;
            if (ss >> x >> y >> z) {
                cfg.extraTargets.push_back({x, y, z});
                double w = 1.0;
                ss >> w;  // опционально; если не прочиталось -- w=1.0
                cfg.extraTargetWeights.push_back(w);
            }
            continue;
        }
        // PATHPT <gid> <x> <y> <z> [weight]  -- точка плотного опорного пути
        // (composite_next). gid тэгирует окно. weight по умолчанию 1.0.
        if (line.rfind("PATHPT ", 0) == 0) {
            std::stringstream ss(line.substr(7));
            int gid; double x, y, z;
            if (ss >> gid >> x >> y >> z) {
                cfg.pathPoints.push_back({x, y, z});
                cfg.pathPointGids.push_back(gid);
                double w = 1.0;
                ss >> w;
                cfg.pathPointWeights.push_back(w);
            }
            continue;
        }
        // GROUP;name;step_abs;step_pct;target;file;value_col;n_matchers;col=val;...
        if (line.rfind("GROUP;", 0) == 0) {
            auto parts = splitDelim(line, ';');
            if (parts.size() < 8 || parts[0] != "GROUP") continue;
            GroupConfig g;
            g.name        = parts[1];
            g.stepAbs     = std::stod(parts[2]);
            g.stepPercent = std::stod(parts[3]);
            g.target      = parts[4];
            g.file        = parts[5];
            g.valueCol    = parts[6];
            int n_match   = std::stoi(parts[7]);
            for (int i = 0; i < n_match && (int)parts.size() > 8 + i; ++i) {
                std::string ms = parts[8 + i];
                auto eq = ms.find('=');
                if (eq == std::string::npos) continue;
                Matcher m;
                m.column = ms.substr(0, eq);
                std::string vals = ms.substr(eq + 1);
                std::stringstream vs(vals);
                std::string v;
                while (std::getline(vs, v, '|')) if (!v.empty()) m.values.push_back(v);
                g.matchers.push_back(std::move(m));
            }
            cfg.groups.push_back(std::move(g));
        }
    }
    return cfg;
}

// ============================================================================
// Метаданные элемента лэттиса (параллельно AcceleratorConfig::elements).
// ============================================================================
struct ElementMeta {
    std::string sp_type;
    int sp_element_id;
    std::string device_type_id;
    int global_id;
    double nominalArg;  // для конверсии arg_val -> currentScale
    double hwCoeff;     // device-level correct_coeff из hardware_library
    double spCoeff;     // per-instance correct_coeff из sp_elements (1.0 default)
                        // эффективное el.correctCoeff = hwCoeff * spCoeff
};

// ============================================================================
// Карта sp_elements: ключ (sp_type, sp_element_id) -> (device_type_id, arg_val).
// Дополнительно храним индекс строки в исходном CSV, чтобы можно было обновлять
// при финальном сохранении.
// ============================================================================
struct SpElementEntry {
    std::string device_type_id;
    double arg_val;
    double correct_coeff;   // per-instance множитель поля (по умолчанию 1.0)
    int row_index;          // индекс в исходном CSV
};

// ============================================================================
// Параметры группы для рантайма (precomputed).
// ============================================================================
struct RuntimeGroup {
    GroupConfig cfg;
    // Индексы MagneticElement в lattice.elements, на которые влияет группа.
    std::vector<int> elementIndices;
    // Текущее значение параметра (одно на всю группу -- единая величина).
    double currentValue = 0.0;
    // Для COEFF/correct_coeff: исходное значение (для DP_* мы НЕ модифицируем
    // sp_elements при сохранении -- это hardware_library).
    // Для arg_val: ключи (sp_type, sp_element_id) для финальной записи обратно.
    // Эти данные понадобятся при сохранении best-state.
    std::vector<std::pair<std::string, int>> spKeys;     // для CURRENT-групп
    std::vector<std::string> hwDeviceTypes;              // для COEFF-групп
};

// ============================================================================
// Helpers для поиска в lattice + sp_elements.
// ============================================================================
// Проверка, удовлетворяет ли (sp_type, sp_element_id, device_type_id) набору матчеров.
static bool matchesAllMatchers(const std::string& sp_type, int sp_eid, const std::string& dev,
                               const std::vector<Matcher>& matchers) {
    for (const auto& m : matchers) {
        std::string actual;
        if      (m.column == "sp_type")        actual = sp_type;
        else if (m.column == "sp_element_id")  actual = std::to_string(sp_eid);
        else if (m.column == "device_type_id") actual = dev;
        else                                   return false; // неизвестная колонка
        bool ok = false;
        for (const auto& v : m.values) if (v == actual) { ok = true; break; }
        if (!ok) return false;
    }
    return true;
}

// ============================================================================
// ЕДИНЫЙ РАСПРЕДЕЛИТЕЛЬ БЮДЖЕТА СИМУЛЯЦИЙ
// ----------------------------------------------------------------------------
// Оба режима (оконный и глобальный line-search) и фаза дрейфа используют ОДИН
// алгоритм: бюджет оценок делится на единицы (группы), и из доли на единицу
// получается число coarse-to-fine раундов на группу. K точек/раунд, потолок
// maxRounds. budget=общий бюджет на данном этапе, units=число групп.
// ============================================================================
static int budgetRounds(int budget, int units, int K, int maxRounds, int minRounds = 1) {
    int r  = std::max(1, budget / std::max(1, units * std::max(1, K)));
    int hi = std::max(1, maxRounds);
    int lo = std::min(std::max(1, minRounds), hi);   // нижний порог не превышает потолок
    // Гарантируем lo раундов на элемент даже при дефиците бюджета (floor),
    // и не больше hi (cap). Бюджетная оценка r определяет точку внутри [lo,hi].
    return std::max(lo, std::min(r, hi));
}

// ============================================================================
// MAIN
// ============================================================================
int main(int argc, char* argv[]) {
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    #endif
    // Безбуферный stdout: важно для отслеживания прогресса в long-running прогонах
    // через background-tasks (пайпы block-буферизуются по умолчанию).
    std::cout.setf(std::ios::unitbuf);

    if (argc < 2) {
        std::cout << "Usage: SynchrotronOptimizer.exe <config.optimconfig>"
                     " [--method=global|windowed|gradient] [--loss=composite|orbit|envelope]\n";
        std::cout << "  Сначала запустите prepare_optim_config.py для конвертации JSON -> optimconfig.\n";
        std::cout << "  --method=global    (default) — ГЛОБАЛЬНЫЙ режим (вся конфигурация сразу);\n";
        std::cout << "                                  метод поиска: global_search_mode=line_search|monte_carlo\n";
        std::cout << "  --method=windowed             — ОКОННЫЙ режим (скользящее окно по магнитам)\n";
        std::cout << "  --method=gradient             — градиентный спуск (FD + Adam + line search)\n";
        std::cout << "  (--mode= — синоним --method=; алиасы: mc/montecarlo->global, per_element/pe->windowed)\n";
        std::cout << "  --budget=N                    — единый бюджет симуляций (override total_iters и n_iter)\n";
        std::cout << "  --ls-points=K                 — line_search: K проб на раунд (override; сим/параметр = 1+rounds*K)\n";
        std::cout << "  --ls-rounds=R                 — line_search: потолок раундов дозума на параметр (override)\n";
        std::cout << "  --search=line_search|monte_carlo|svd — метод поиска в global-режиме (override JSON)\n";
        std::cout << "  --kick-every=N --kick-scale=A — monte_carlo: частота/масштаб киков (override; 0=без киков)\n";
        std::cout << "  --loss=composite   (default) — position + direction + envelope\n";
        std::cout << "  --loss=orbit                  — только position (sum d² от опорной орбиты),\n";
        std::cout << "                                  для грубой настройки\n";
        std::cout << "  --loss=envelope               — только envelope (рост бетатронных амплитуд),\n";
        std::cout << "                                  для подавления колебаний при готовой орбите\n";
        return 1;
    }
    std::string ws = findWorkspace(argv[0]);
    std::cout << "Workspace: " << ws << "\n";

    StepTimer totalTimer;

    // ---- 1. Load flat config ----
    std::string cfgPath = argv[1];
    OptimConfig cfg = loadFlatConfig(cfgPath);
    if (cfg.expId == 0 || cfg.groups.empty()) {
        std::cerr << "ERROR: empty config or no enabled groups in " << cfgPath << "\n";
        return 1;
    }
    // CLI-флаги переопределяют JSON.
    //   --method=mc|grad
    //   --loss=composite|orbit
    for (int ai = 2; ai < argc; ++ai) {
        std::string a = argv[ai];
        if (a.rfind("--method=", 0) == 0 || a.rfind("--mode=", 0) == 0) {
            std::string v = a.substr(a.find('=') + 1);
            std::string cm = canonicalMode(v);
            if (cm == "global" || cm == "windowed" || cm == "lens" || cm == "gradient") {
                cfg.method = cm;
                // --method=mc/montecarlo/random исторически = Monte Carlo: ставим и метод поиска,
                // иначе остаётся line_search по умолчанию (частая путаница). --search= перебивает.
                if (canonicalSearchMethod(v) == "monte_carlo") cfg.globalSearchMode = "monte_carlo";
            }
            else { std::cerr << "ERROR: --method/--mode must be global|windowed|lens|gradient "
                                "(алиасы: mc/montecarlo->global, per_element/pe->windowed, "
                                "grad->gradient), got '" << v << "'\n"; return 1; }
        }
        else if (a.rfind("--loss=", 0) == 0) {
            std::string v = a.substr(7);
            if      (v == "orbit" || v == "rough")        cfg.lossMode = "orbit";
            else if (v == "composite" || v == "full")     cfg.lossMode = "composite";
            else if (v == "envelope" || v == "betatron")  cfg.lossMode = "envelope";
            else { std::cerr << "ERROR: --loss must be orbit|composite|envelope, got '" << v << "'\n"; return 1; }
        }
        // Per-element быстрые override'ы из CLI (без правки JSON):
        else if (a.rfind("--W=", 0) == 0)         cfg.perElement.windowSize     = std::stoi(a.substr(4));
        else if (a.rfind("--stride=", 0) == 0)    cfg.perElement.windowStride   = std::stoi(a.substr(9));
        else if (a.rfind("--iters=", 0) == 0)     cfg.perElement.itersPerWindow = std::stoi(a.substr(8));
        else if (a.rfind("--passes=", 0) == 0)    cfg.perElement.passes         = std::stoi(a.substr(9));
        else if (a.rfind("--n=", 0) == 0)         cfg.nIter                     = std::stoi(a.substr(4));
        // --budget= : ЕДИНЫЙ бюджет симуляций. Переопределяет total_iters (бюджет
        // line-search в оконном и глобальном режимах) И n_iter (для метода
        // monte_carlo). Алиас --total_iters=.
        else if (a.rfind("--budget=", 0) == 0 || a.rfind("--total_iters=", 0) == 0) {
            int b = std::stoi(a.substr(a.find('=') + 1));
            cfg.perElement.totalIters = b;
            cfg.nIter = b;
        }
        // line-search: K проб на раунд и потолок раундов на параметр (override JSON).
        // Симуляций на 1 параметр за свип = 1 (incumbent) + rounds × K.
        else if (a.rfind("--ls-points=", 0) == 0)  cfg.perElement.lineSearchPoints    = std::stoi(a.substr(12));
        else if (a.rfind("--ls-rounds=", 0) == 0)  cfg.perElement.lineSearchMaxRounds = std::stoi(a.substr(12));
        // --search= : метод поиска в global-режиме (line_search | monte_carlo | svd), override JSON.
        else if (a.rfind("--search=", 0) == 0)     cfg.globalSearchMode = canonicalSearchMethod(a.substr(9));
        // monte_carlo: частота и масштаб киков (override JSON). 0 = без киков.
        else if (a.rfind("--kick-every=", 0) == 0) cfg.kickEvery = std::stoi(a.substr(13));
        else if (a.rfind("--kick-scale=", 0) == 0) cfg.kickScale = std::stod(a.substr(13));
    }

    // Применение режима loss.
    // - "orbit":    отключаем direction и envelope (только position).
    // - "envelope": отключаем position и direction (только envelope).
    // - "composite": все три компоненты, как задано в JSON.
    // Оригинальные веса сохраняем для информативного лога.
    double origDirWeight = cfg.dirWeight;
    double origEnvWeight = cfg.envWeight;
    if (cfg.lossMode == "orbit") {
        cfg.dirWeight = 0.0;
        cfg.envWeight = 0.0;
    } else if (cfg.lossMode == "envelope") {
        cfg.usePositionLoss = false;
        cfg.dirWeight = 0.0;
        if (cfg.envWeight <= 0.0) {
            std::cerr << "ERROR: --loss=envelope требует envelope_loss.weight > 0 в JSON\n";
            return 1;
        }
    }
    std::cout << "Config: " << cfgPath << "\n";
    std::cout << "  exp_id=" << cfg.expId
              << "  mode=" << cfg.method
              << "  loss=" << cfg.lossMode;
    if (cfg.lossMode == "orbit" && (origDirWeight > 0 || origEnvWeight > 0)) {
        std::cout << " [direction(w=" << origDirWeight
                  << ") и envelope(w=" << origEnvWeight << ") ОТКЛЮЧЕНЫ]";
    } else if (cfg.lossMode == "envelope") {
        std::cout << " [position и direction(w=" << origDirWeight << ") ОТКЛЮЧЕНЫ,"
                  << " envelope(w=" << origEnvWeight << ") активен]";
    }
    std::cout << "  groups=" << cfg.groups.size();
    if (cfg.method == "global") {
        std::cout << "  n_iter=" << cfg.nIter
                  << "  kick_every=" << cfg.kickEvery
                  << "  kick_scale=" << cfg.kickScale;
    } else {
        std::cout << "  n_iter=" << cfg.gradient.nIter
                  << "  lr=" << cfg.gradient.learningRate
                  << "  fd=" << cfg.gradient.fdMethod
                  << "  ε_abs=" << cfg.gradient.fdStepAbs
                  << "  ε_pct=" << cfg.gradient.fdStepPercent;
    }
    std::cout << "\n";

    // ---- 2. Load experiment ----
    auto expCsv = readCsvTable(ws + "tables/experiments.csv");
    int latId = -1, grpId = -1, maxTurns = 0;
    int trajStride = 1;
    double dt = 0.0;
    FieldVector fp1 = {5.0, -1.0, -1.0};
    FieldVector fp2 = {5.0,  1.0, -1.0};
    FieldVector fp3 = {5.0,  0.0,  1.0};
    FieldVector fpDir = {0.0, 0.0, 0.0};   // направление пересечения (опц., колонки 15-17)
    for (auto& r : expCsv.rows) {
        if ((int)r.size() < 5) continue;
        if (std::stoi(r[0]) == cfg.expId) {
            latId    = std::stoi(r[1]);
            grpId    = std::stoi(r[2]);
            dt       = std::stod(r[3]);
            maxTurns = std::stoi(r[4]);
            // Колонка 5 -- traj_stride.
            if (r.size() >= 15) {
                try { trajStride = std::max(1, std::stoi(r[5])); } catch (...) {}
                fp1 = {std::stod(r[6]),  std::stod(r[7]),  std::stod(r[8])};
                fp2 = {std::stod(r[9]),  std::stod(r[10]), std::stod(r[11])};
                fp3 = {std::stod(r[12]), std::stod(r[13]), std::stod(r[14])};
            }
            if (r.size() >= 18) {   // направление пересечения финиш-плоскости
                fpDir = {std::stod(r[15]), std::stod(r[16]), std::stod(r[17])};
            }
            break;
        }
    }
    if (latId < 0) {
        std::cerr << "ERROR: experiment " << cfg.expId << " not found\n";
        return 1;
    }
    std::cout << "  experiment: lat=" << latId << " grp=" << grpId
              << " dt=" << dt << " skipSteps=" << maxTurns << "\n";

    // ---- 3. Load all CSVs into memory ----
    auto spCsv     = readCsvTable(ws + "tables/sp_elements.csv");
    auto layoutCsv = readCsvTable(ws + "tables/lattice_layout.csv");
    auto hwCsv     = readCsvTable(ws + "tables/hardware_library.csv");

    // sp_elements lookup map keyed by (sp_type, sp_element_id)
    int sp_idx_type  = spCsv.col("sp_type");
    int sp_idx_eid   = spCsv.col("sp_element_id");
    int sp_idx_dev   = spCsv.col("device_type_id");
    int sp_idx_val   = spCsv.col("arg_val");
    int sp_idx_coeff = spCsv.col("correct_coeff");   // -1 если колонки нет
    std::map<std::pair<std::string, int>, SpElementEntry> spMap;
    for (int i = 0; i < (int)spCsv.rows.size(); ++i) {
        auto& r = spCsv.rows[i];
        if ((int)r.size() <= sp_idx_val) continue;
        try {
            SpElementEntry e;
            e.device_type_id = r[sp_idx_dev];
            // Пустой arg_val == "-" == 0.0: постоянные магниты (DP_*_PM, Active=False)
            // не имеют тока, поле берётся из correct_coeff. БАГ (исправлено): раньше
            // std::stod("") бросал исключение -> вся строка отбрасывалась из spMap, и ВСЕ
            // PM-диполи выпадали из решётки оптимизатора (beam без поворотных!), хотя в
            // трекере (lattice_configs.csv) они есть. Теперь оптимизатор согласован с трекером.
            e.arg_val        = (r[sp_idx_val] == "-" || r[sp_idx_val].empty())
                                   ? 0.0 : std::stod(r[sp_idx_val]);
            // Per-instance correct_coeff (опционально, default 1.0)
            e.correct_coeff  = 1.0;
            if (sp_idx_coeff >= 0 && (int)r.size() > sp_idx_coeff &&
                !r[sp_idx_coeff].empty() && r[sp_idx_coeff] != "-") {
                try { e.correct_coeff = std::stod(r[sp_idx_coeff]); } catch (...) {}
            }
            e.row_index      = i;
            std::string spt  = r[sp_idx_type];
            int seid         = std::stoi(r[sp_idx_eid]);
            spMap[{spt, seid}] = e;
        } catch (...) {}
    }
    std::cout << "  sp_elements: " << spMap.size() << " entries\n";

    // ---- 4. Build lattice (assemble equivalent, in-memory) ----
    StepTimer initTimer;
    AcceleratorConfig lattice;
    std::vector<ElementMeta> meta;
    std::map<std::string, std::shared_ptr<FieldMap>> mapCache;

    int lay_gid  = layoutCsv.col("global_id");
    int lay_spn  = layoutCsv.col("sp_number");
    int lay_typ  = layoutCsv.col("sp_type");
    int lay_eid  = layoutCsv.col("sp_element_id");
    int lay_x    = layoutCsv.col("x");
    int lay_y    = layoutCsv.col("y");
    int lay_z    = layoutCsv.col("z");
    int lay_yaw  = layoutCsv.col("yaw");
    int lay_pit  = layoutCsv.col("pitch");
    int lay_roll = layoutCsv.col("roll");

    int hw_id    = hwCsv.col("device_type_id");
    int hw_act   = hwCsv.col("Active");
    int hw_bin   = hwCsv.col("map_bin_path");
    int hw_fd    = hwCsv.col("field_direction");
    int hw_cc    = hwCsv.col("correct_coeff");
    int hw_ft    = hwCsv.col("field_type");      // -1 если колонка отсутствует
    int hw_sym   = hwCsv.col("use_symmetry");    // -1 если колонка отсутствует

    auto findHwRow = [&](const std::string& devType) -> int {
        for (int i = 0; i < (int)hwCsv.rows.size(); ++i)
            if ((int)hwCsv.rows[i].size() > hw_id && hwCsv.rows[i][hw_id] == devType) return i;
        return -1;
    };

    int n_skipped_by_gid = 0;
    for (auto& r : layoutCsv.rows) {
        if ((int)r.size() <= lay_roll) continue;
        // Если задан latticeMaxGid, отсекаем элементы с global_id больше.
        // Это нужно для экспериментов, где трекинг идёт через подсет lattice
        // (например, 3 SP + свободный пролёт за ними).
        if (cfg.latticeMaxGid > 0) {
            int gid;
            try { gid = std::stoi(r[lay_gid]); } catch (...) { continue; }
            if (gid > cfg.latticeMaxGid) { n_skipped_by_gid++; continue; }
        }
        std::string spt = r[lay_typ];
        int seid;
        try { seid = std::stoi(r[lay_eid]); } catch (...) { continue; }
        auto it = spMap.find({spt, seid});
        if (it == spMap.end()) continue;
        const std::string& devType = it->second.device_type_id;
        double argVal = it->second.arg_val;

        int hi = findHwRow(devType);
        if (hi < 0) continue;
        const auto& h = hwCsv.rows[hi];
        bool isActive = (h[hw_act] == "True" || h[hw_act] == "true" || h[hw_act] == "1");
        std::string binName = h[hw_bin];
        double fieldDir = std::stod(h[hw_fd]);
        double corrCoeff = std::stod(h[hw_cc]);

        // field_type: "B"/"M"/"MAGNETIC" → MAGNETIC; "E"/"ELECTRIC" → ELECTRIC.
        // Колонка отсутствует или значение нестандартное → MAGNETIC.
        FieldType fieldKind = FieldType::MAGNETIC;
        if (hw_ft >= 0 && (int)h.size() > hw_ft) {
            std::string ft = h[hw_ft];
            for (auto& c : ft) c = (char)std::toupper((unsigned char)c);
            if (ft == "E" || ft == "ELECTRIC") fieldKind = FieldType::ELECTRIC;
        }
        // use_symmetry: "1"/"True" → true, "0"/"False" → false, иначе → true.
        bool useSym = true;
        if (hw_sym >= 0 && (int)h.size() > hw_sym && !h[hw_sym].empty()) {
            useSym = (h[hw_sym] == "1" || h[hw_sym] == "True" || h[hw_sym] == "true");
        }

        if (mapCache.find(devType) == mapCache.end()) {
            auto m = std::make_shared<FieldMap>();
            if (!m->loadFromBinary(ws + "fields_bin/" + binName)) continue;
            mapCache[devType] = m;
        }
        MagneticElement el;
        el.map            = mapCache[devType];
        el.type           = fieldKind;
        el.useSymmetry    = useSym;
        el.isActive       = isActive;
        el.fieldDirection = fieldDir;
        // Эффективный коэффициент = hw × per-instance из sp_elements.
        el.correctCoeff   = corrCoeff * it->second.correct_coeff;
        el.setOrientation(std::stod(r[lay_x]),  std::stod(r[lay_y]),  std::stod(r[lay_z]),
                          std::stod(r[lay_yaw]), std::stod(r[lay_pit]), std::stod(r[lay_roll]));
        if (isActive) el.currentScale = argVal / el.map->getNominalArg();
        else          el.currentScale = 1.0;

        ElementMeta em;
        em.sp_type        = spt;
        em.sp_element_id  = seid;
        em.device_type_id = devType;
        em.global_id      = std::stoi(r[lay_gid]);
        em.nominalArg     = el.map->getNominalArg();
        em.hwCoeff        = corrCoeff;
        em.spCoeff        = it->second.correct_coeff;
        meta.push_back(em);
        lattice.addElement(el);
    }
    std::cout << "  lattice elements: " << lattice.size()
              << " (" << mapCache.size() << " unique types)"
              << "  init=" << initTimer.ms() << " ms\n";
    if (cfg.latticeMaxGid > 0) {
        std::cout << "  LATTICE_MAX_GID=" << cfg.latticeMaxGid
                  << " -> skipped " << n_skipped_by_gid << " elements past gid limit\n";
    }

    // Spatial index для O(1) поиска магнитов в окрестности позиции.
    {
        StepTimer t_idx;
        int totalRegs = lattice.buildSpatialIndex();
        int nbx, nby, nNonEmpty, maxLen; double bsize, avgLen;
        lattice.getGridStats(nbx, nby, bsize, nNonEmpty, maxLen, avgLen);
        std::cout << "  spatial index: " << t_idx.ms() << " ms, "
                  << nbx << "x" << nby << " buckets (size=" << bsize << " m), "
                  << nNonEmpty << " non-empty, avg=" << avgLen
                  << "/bucket, max=" << maxLen
                  << " (" << totalRegs << " regs)\n";
    }

    // Дизайн-полилиния = упорядоченный список позиций элементов lattice
    // + одна точка экстраполяции на каждый конец (prelude/postlude по 10/20 м).
    // Используется в envelope_loss для измерения отклонения частиц от опорной
    // орбиты (включая P1, которая в нашей симуляции тоже может уплывать).
    std::vector<FieldVector> designPolyline;
    if (cfg.envWeight > 0.0 && !lattice.elementsRef().empty()) {
        const auto& els = lattice.elementsRef();
        std::vector<FieldVector> body;
        std::string source = "elements";

        // Приоритет 1: готовый CSV из prepare_optim_config.py (истинная
        // проектная орбита, собранная из .lte).
        if (!cfg.designOrbitCsv.empty()) {
            body = loadDesignOrbitCsv(cfg.designOrbitCsv);
            if (body.size() >= 2) {
                source = "lte:" + cfg.designOrbitCsv;
            } else {
                std::cerr << "  WARN: DESIGN_ORBIT_CSV не загружается, "
                             "откатываюсь к построению через центры элементов\n";
                body.clear();
            }
        }

        // Fallback: старое поведение -- центры элементов (опционально с дугами).
        if (body.empty()) {
            std::vector<FieldVector> centers;
            std::vector<double>      yaws;
            centers.reserve(els.size());
            yaws.reserve(els.size());
            for (const auto& el : els) {
                centers.push_back({el.transform.tx, el.transform.ty, el.transform.tz});
                yaws.push_back(std::atan2(el.transform.m[1][0], el.transform.m[0][0]));
            }
            if (cfg.envPolylineMode == "arcs") {
                body = refinePolylineWithArcs(centers, yaws, cfg.envArcSubdivide);
                source = "centers:arcs";
            } else {
                body = centers;
                source = "centers:linear";
            }
        }

        // prelude/postlude -- линейные продолжения по first/last tangent.
        // Касательные берём из first/last элементов lattice (валидно и при
        // CSV-источнике: orbit стартует с (0,0,0), элементы тоже).
        const auto& first = els.front();
        const auto& last  = els.back();
        FieldVector t0 = {first.transform.m[0][0],
                          first.transform.m[1][0],
                          first.transform.m[2][0]};
        FieldVector tN = {last.transform.m[0][0],
                          last.transform.m[1][0],
                          last.transform.m[2][0]};
        const FieldVector& first_pos = body.front();
        const FieldVector& last_pos  = body.back();
        designPolyline.reserve(body.size() + 2);
        designPolyline.push_back({first_pos.x - 10.0*t0.x,
                                  first_pos.y - 10.0*t0.y,
                                  first_pos.z - 10.0*t0.z});
        for (const auto& p : body) designPolyline.push_back(p);
        designPolyline.push_back({last_pos.x + 20.0*tN.x,
                                  last_pos.y + 20.0*tN.y,
                                  last_pos.z + 20.0*tN.z});
        std::cout << "  design polyline: " << designPolyline.size()
                  << " pts (source=" << source << ", +2 extensions)\n";
    }

    // ---- 5. Beam ----
    auto pgCsv = readCsvTable(ws + "tables/particle_groups.csv");
    std::vector<ParticleState> beam, beam0;
    std::vector<double> masses, charges;
    for (auto& r : pgCsv.rows) {
        if ((int)r.size() < 10) continue;
        if (std::stoi(r[0]) != grpId) continue;
        masses.push_back(std::stod(r[2]));
        charges.push_back(std::stod(r[3]));
        ParticleState s;
        s.pos = {std::stod(r[4]), std::stod(r[5]), std::stod(r[6])};
        FieldVector v0 = {std::stod(r[7]), std::stod(r[8]), std::stod(r[9])};
        s.mom = Integrator::velocityToMomentum(v0, masses.back());
        beam.push_back(s);
        beam0.push_back(s);
    }
    // Ограничение числа частиц (фаза поворотных — только центральная P1).
    if (cfg.maxParticles > 0 && (int)beam.size() > cfg.maxParticles) {
        int keep = cfg.maxParticles;
        beam.resize(keep); beam0.resize(keep);
        masses.resize(keep); charges.resize(keep);
        std::cout << "  particles: " << beam.size()
                  << " (ограничено MAX_PARTICLES=" << cfg.maxParticles << ")\n";
    } else {
        std::cout << "  particles: " << beam.size() << "\n";
    }

    // ---- 6. Magnet positions for loss ----
    // ИСТОЧНИК КООРДИНАТ:
    //   Используем середины (entry+exit)/2 из magnet_endpoints_lat{N}.csv --
    //   это точки на ИСТИННОЙ эталонной орбите (собранной из .lte). Это
    //   важно: координаты x,y в lattice_layout.csv -- центры ЖЕЛЕЗА магнитов,
    //   а они в дипольных секциях смещены от орбиты на sagitta ≈ до 17 мм.
    //   До этой правки position-loss мерял до центров железа, что давало
    //   ложно-завышенный mean при правильно настроенной орбите.
    //   Fallback на lattice_layout остаётся для магнитов, для которых
    //   endpoint не найден (например, при отсутствии CSV).
    std::unordered_map<int, MagnetEndpoint> endpointsByGidGlobal;
    if (!cfg.magnetEndpointsCsv.empty()) {
        auto eps = loadMagnetEndpointsCsv(cfg.magnetEndpointsCsv);
        for (const auto& e : eps) endpointsByGidGlobal[e.global_id] = e;
        std::cout << "  endpoints loaded: " << endpointsByGidGlobal.size()
                  << " (для PE start_pos и envelope_loss)\n";
    }
    // Position-loss reference: iron (default, центры железа) или orbit (середина дуги).
    const bool use_orbit_for_pos = (cfg.positionReference == "orbit");
    std::cout << "  position-loss reference: "
              << (use_orbit_for_pos ? "orbit (центр эталонной орбиты)"
                                    : "iron (центры железа из lattice_layout)")
              << "\n";

    std::set<int> spFilter(cfg.lossSpNumbers.begin(), cfg.lossSpNumbers.end());
    std::set<std::string> devExclude(cfg.lossExcludeDevices.begin(),
                                     cfg.lossExcludeDevices.end());
    std::vector<FieldVector> magnetPos;
    // Параллельный набор АЛЬТЕРНАТИВНОГО референса (для диагностики dual-loss):
    // если оптимизируем по железу -> здесь .lte-орбита, и наоборот. Тот же
    // порядок индексов, что и magnetPos.
    std::vector<FieldVector> magnetPosAlt;
    std::vector<int>         magnetSp;
    // Параллельный массив global_id (для фильтрации в PE composite режиме).
    // -1 для extra TARGET'ов (не имеют lattice-gid). Извне магниты идут с gid > 0.
    std::vector<int>         magnetGid;
    // Параллельный массив весов: 1.0 для обычных магнитов из layout,
    // настраиваемое значение для extra TARGET'ов (см. строки TARGET в optimconfig).
    std::vector<double>      magnetWeight;
    std::set<std::tuple<double, double, double>> seenPos;
    int n_excluded_by_device = 0;
    for (auto& r : layoutCsv.rows) {
        if ((int)r.size() <= lay_z) continue;
        // Тот же gid-фильтр, что и при сборке lattice -- чтобы в loss не попали
        // магниты-цели, которые не "существуют" в текущей урезанной решётке.
        if (cfg.latticeMaxGid > 0) {
            int gid;
            try { gid = std::stoi(r[lay_gid]); } catch (...) { continue; }
            if (gid > cfg.latticeMaxGid) continue;
        }
        int spn;
        try { spn = std::stoi(r[lay_spn]); } catch (...) { continue; }
        if (!spFilter.empty() && !spFilter.count(spn)) continue;
        // Фильтр по device_type_id: lookup через spMap по (sp_type, sp_element_id).
        if (!devExclude.empty()) {
            std::string spt = r[lay_typ];
            int seid;
            try { seid = std::stoi(r[lay_eid]); } catch (...) { continue; }
            auto it = spMap.find({spt, seid});
            if (it != spMap.end() && devExclude.count(it->second.device_type_id)) {
                n_excluded_by_device++;
                continue;
            }
        }
        // Координаты для loss -- считаем ОБА референса:
        //   iron  -> центр железа из lattice_layout.
        //   orbit -> середина (entry+exit)/2 из magnet_endpoints (.lte-орбита).
        // primary = по configu (use_orbit_for_pos), alt = другой (для диагностики).
        int gid_for_pos = -1;
        try { gid_for_pos = std::stoi(r[lay_gid]); } catch (...) {}
        double x_iron = std::stod(r[lay_x]);
        double y_iron = std::stod(r[lay_y]);
        double z_iron = std::stod(r[lay_z]);
        double x_orb = x_iron, y_orb = y_iron, z_orb = z_iron;  // fallback = iron
        if (gid_for_pos > 0) {
            auto eit = endpointsByGidGlobal.find(gid_for_pos);
            if (eit != endpointsByGidGlobal.end()) {
                const auto& ep = eit->second;
                x_orb = 0.5 * (ep.x_in + ep.x_out);
                y_orb = 0.5 * (ep.y_in + ep.y_out);
                z_orb = 0.0;
            }
        }
        double x, y, z, xa, ya, za;
        if (use_orbit_for_pos) { x=x_orb;  y=y_orb;  z=z_orb;
                                 xa=x_iron; ya=y_iron; za=z_iron; }
        else                   { x=x_iron; y=y_iron; z=z_iron;
                                 xa=x_orb;  ya=y_orb;  za=z_orb; }
        auto key = std::make_tuple(std::round(x*1e6)/1e6,
                                   std::round(y*1e6)/1e6,
                                   std::round(z*1e6)/1e6);
        if (seenPos.count(key)) continue;
        seenPos.insert(key);
        magnetPos.push_back({x, y, z});
        magnetPosAlt.push_back({xa, ya, za});
        magnetSp.push_back(spn);
        magnetGid.push_back(gid_for_pos);
        magnetWeight.push_back(1.0);  // обычные магниты получают вес 1.0
    }
    std::cout << "  magnets (loss, sp_numbers=";
    for (int n : cfg.lossSpNumbers) std::cout << n << ",";
    std::cout << "): " << magnetPos.size() << "\n";
    if (!devExclude.empty()) {
        std::cout << "  excluded by device_type (";
        for (const auto& d : cfg.lossExcludeDevices) std::cout << d << ",";
        std::cout << "): " << n_excluded_by_device << " rows\n";
    }

    // Дополнительные explicit-цели (например, центр первого магнита SP5).
    // Для extra targets дедупликация ОТКЛЮЧЕНА -- допускаем несколько TARGET'ов
    // в одной точке, это эквивалентно увеличению weight (но менее удобно).
    int n_extra = 0;
    double sum_extra_weight = 0.0;
    for (size_t ti = 0; ti < cfg.extraTargets.size(); ++ti) {
        const auto& t = cfg.extraTargets[ti];
        double w = (ti < cfg.extraTargetWeights.size()) ? cfg.extraTargetWeights[ti] : 1.0;
        magnetPos.push_back(t);
        magnetPosAlt.push_back(t);  // extra-цель одинакова в обоих референсах
        magnetSp.push_back(-1);  // -1 = explicit target (не из sp_number-фильтра)
        magnetGid.push_back(-1); // -1 = extra TARGET, всегда включается в окно
        magnetWeight.push_back(w);
        n_extra++;
        sum_extra_weight += w;
    }
    if (n_extra > 0) {
        std::cout << "  extra targets added: " << n_extra
                  << " (sum weight=" << sum_extra_weight
                  << ", total magnets for loss: " << magnetPos.size() << ")\n";
    }

    // ---- 7. Precompute group -> element indices ----
    std::vector<RuntimeGroup> rgroups;
    for (auto& g : cfg.groups) {
        RuntimeGroup rg;
        rg.cfg = g;
        // Найти MagneticElement, удовлетворяющие всем матчерам.
        for (int ei = 0; ei < (int)meta.size(); ++ei) {
            if (matchesAllMatchers(meta[ei].sp_type, meta[ei].sp_element_id,
                                   meta[ei].device_type_id, g.matchers)) {
                rg.elementIndices.push_back(ei);
            }
        }
        // Текущее значение параметра -- считываем из meta[].
        if (g.target == "COEFF_HW") {
            // hardware_library.correct_coeff: одно значение на device_type, общее.
            if (!rg.elementIndices.empty())
                rg.currentValue = meta[rg.elementIndices.front()].hwCoeff;
            for (const auto& m : g.matchers)
                if (m.column == "device_type_id")
                    for (const auto& v : m.values) rg.hwDeviceTypes.push_back(v);
        } else if (g.target == "COEFF_INST") {
            // sp_elements.correct_coeff: per-instance, идентифицируется по
            // (sp_type, sp_element_id).
            if (!rg.elementIndices.empty())
                rg.currentValue = meta[rg.elementIndices.front()].spCoeff;
            std::set<std::pair<std::string, int>> keys;
            for (int ei : rg.elementIndices)
                keys.insert({meta[ei].sp_type, meta[ei].sp_element_id});
            for (auto& k : keys) rg.spKeys.push_back(k);
        } else {  // CURRENT
            if (!rg.elementIndices.empty()) {
                const auto& el0 = lattice.elementsRef()[rg.elementIndices.front()];
                rg.currentValue = el0.currentScale * meta[rg.elementIndices.front()].nominalArg;
            }
            std::set<std::pair<std::string, int>> keys;
            for (int ei : rg.elementIndices)
                keys.insert({meta[ei].sp_type, meta[ei].sp_element_id});
            for (auto& k : keys) rg.spKeys.push_back(k);
        }
        if (rg.elementIndices.empty()) {
            std::cout << "  WARN: group '" << g.name << "' matched 0 elements (skipped)\n";
            continue;
        }
        rgroups.push_back(std::move(rg));
    }
    std::cout << "  active groups: " << rgroups.size() << "\n";

    // ---- 7b. Пометка групп-линз для режима lagged-lens ----
    // groupLagged[gi] == 1 -> группа фокусирующая (device_type в lensLaggedDeviceTypes),
    // активируется с задержкой lensLagWindows окон. Иначе (поворотные) -- сразу.
    std::vector<char> groupLagged(rgroups.size(), 0);
    if (!cfg.perElement.lensLaggedDeviceTypes.empty()) {
        int n_lagged = 0;
        for (size_t gi = 0; gi < rgroups.size(); ++gi) {
            const auto& ei = rgroups[gi].elementIndices;
            if (!ei.empty() &&
                cfg.perElement.lensLaggedDeviceTypes.count(meta[ei.front()].device_type_id)) {
                groupLagged[gi] = 1;
                ++n_lagged;
            }
        }
        if (cfg.perElement.lensLagWindows > 0)
            std::cout << "  lagged-lens: " << n_lagged << " групп-линз, лаг "
                      << cfg.perElement.lensLagWindows << " окон\n";
    }

    // ---- 8. Apply value to lattice (mutate currentScale / correctCoeff) ----
    // Эффективный el.correctCoeff = meta[ei].hwCoeff * meta[ei].spCoeff.
    // Перенастраиваем соответствующий слот в meta[] и пересчитываем эффект.
    auto applyGroupValue = [&](RuntimeGroup& rg, double value) {
        rg.currentValue = value;
        auto& els = lattice.elementsMut();
        if (rg.cfg.target == "COEFF_HW") {
            for (int ei : rg.elementIndices) {
                meta[ei].hwCoeff = value;
                els[ei].correctCoeff = meta[ei].hwCoeff * meta[ei].spCoeff;
            }
        } else if (rg.cfg.target == "COEFF_INST") {
            for (int ei : rg.elementIndices) {
                meta[ei].spCoeff = value;
                els[ei].correctCoeff = meta[ei].hwCoeff * meta[ei].spCoeff;
            }
        } else {  // CURRENT
            for (int ei : rg.elementIndices) {
                if (els[ei].isActive) els[ei].currentScale = value / meta[ei].nominalArg;
            }
        }
    };

    // ---- 9. Tracking + loss ----
    FieldVector startRef = beam.empty() ? FieldVector{0,0,0} : beam[0].pos;

    // ИНКРЕМЕНТАЛЬНЫЙ LOSS: вместо того, чтобы складывать всю траекторию в память,
    // а потом строить воксельную сетку и искать ближайшую точку для каждого магнита,
    // мы прямо в цикле трекинга поддерживаем minDist²[i][m] для каждой пары
    // (частица, магнит). На каждом шаге обновляем min — это всего ~8 flops × M магнитов
    // на шаг. На выходе: loss = sum(d²) уже посчитан, kdtree-проход больше не нужен.
    //
    // Профит: (1) нет push_back на каждом шаге → освобождается L3 cache под field maps,
    //          (2) пропускаем воксельную сетку (~150-200 ms / итерацию),
    //          (3) loss считается по ВСЕМ шагам интегрирования (с trajStride=1
    //              старого подхода тоже), но БЕЗ хранения 1.5M точек в памяти.
    //
    // Возвращает true, если все 5 частиц достигли финиша. Заполняет
    // outLoss / outMeanMM / outMaxMM (тот же формат, что был у computeLoss).
    // Расширенная подпись: опциональные указатели outPosComp/outDirComp/outEnvComp
    // возвращают компоненты loss БЕЗ применения базовых весов (cfg.dirWeight,
    // cfg.envWeight). Это нужно для dynamic_loss reweighting в MC: одни и те же
    // компоненты переvзвешиваются разными весами на каждой итерации.
    // outLoss всегда содержит СТАТИЧЕСКИЙ composite (как раньше): pos +
    // cfg.dirWeight·dir + cfg.envWeight·env. Поведение для callers без
    // указателей -- bit-for-bit идентично прежнему.
    auto runTrackingAndLoss = [&](double& outLoss, double& outMeanMM, double& outMaxMM,
                                   double* outPosComp = nullptr,
                                   double* outDirComp = nullptr,
                                   double* outEnvComp = nullptr,
                                   double* outParallelComp = nullptr,
                                   double* outDecayComp = nullptr,
                                   std::vector<double>* outResidual = nullptr) -> bool {
        for (size_t i = 0; i < beam.size(); ++i) beam[i] = beam0[i];
        std::vector<FinishPlane> planes(beam.size());
        for (size_t i = 0; i < beam.size(); ++i)
            planes[i].setup(fp1, fp2, fp3, startRef, maxTurns, fpDir);

        const size_t M = magnetPos.size();
        const double INF = std::numeric_limits<double>::infinity();
        std::vector<std::vector<double>> minD2(beam.size(), std::vector<double>(M, INF));
        std::vector<char> reached(beam.size(), 0);
        // Меж-частичный член: позиции ближайшего подхода (только если включён).
        const bool interOn = (cfg.interparticleWeight > 0.0);
        const bool decayOn = (cfg.envelopeDecayWeight > 0.0);
        const bool needClose = interOn || decayOn;
        std::vector<std::vector<FieldVector>> closePos;
        if (needClose) closePos.assign(beam.size(), std::vector<FieldVector>(M));

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < (int)beam.size(); ++i) {
            auto& md = minD2[i];
            // Стартовая точка -- сразу обновляем minDist²
            {
                const FieldVector& p = beam[i].pos;
                for (size_t m = 0; m < M; ++m) {
                    double dx = p.x - magnetPos[m].x;
                    double dy = p.y - magnetPos[m].y;
                    double dz = p.z - magnetPos[m].z;
                    md[m] = dx*dx + dy*dy + dz*dz;
                    if (needClose) closePos[i][m] = p;
                }
            }
            int steps = 0;
            while (steps < 50000000) {
                Integrator::borisStep(beam[i], masses[i], charges[i], lattice, dt);
                steps++;
                // Обновить min² до всех магнитов в текущей позиции.
                const FieldVector& p = beam[i].pos;
                const FieldVector* mp = magnetPos.data();
                for (size_t m = 0; m < M; ++m) {
                    double dx = p.x - mp[m].x;
                    double dy = p.y - mp[m].y;
                    double dz = p.z - mp[m].z;
                    double d2 = dx*dx + dy*dy + dz*dz;
                    if (d2 < md[m]) { md[m] = d2; if (needClose) closePos[i][m] = p; }
                }
                if (planes[i].hasCrossed(beam[i].pos)) { reached[i] = 1; break; }
            }
        }
        // Все 5 частиц обязаны были достигнуть финиша.
        for (size_t i = 0; i < beam.size(); ++i) if (!reached[i]) return false;

        // Суммируем loss. ВЕС из magnetWeight применяется только к loss-сумме;
        // mean/max остаются "геометрическими" статистиками (без веса) -- они
        // нужны как диагностика средней/максимальной невязки в мм.
        double weighted_sum_d2 = 0.0, weighted_sum_d = 0.0;
        double sum_d = 0.0, max_d = 0.0;
        size_t total = 0;
        // outResidual (для SVD/Гаусс-Ньютона): вектор невязок r_{i,m} = sqrt(w_m)·d(i,m),
        // так что Σ r² = Σ w_m·d² = position loss (squared). Минимизация ||r||² = коррекция
        // позиции пучка к центрам магнитов (классическая коррекция орбиты). При
        // svdAngleWeight>0 в хвост добавляются 2 ПЕРИОДИЧЕСКИЕ записи ЦЕНТРОИДА (см. ниже).
        const bool angleActive = (cfg.svdAngleWeight > 0.0 &&
            (fpDir.x*fpDir.x + fpDir.y*fpDir.y + fpDir.z*fpDir.z) > 1e-24);
        const bool angleRes = (outResidual && angleActive);   // только когда нужен вектор невязок
        if (outResidual) outResidual->assign(beam.size() * M + (angleRes ? 2 : 0), 0.0);
        for (size_t i = 0; i < minD2.size(); ++i) {
            const auto& md = minD2[i];
            for (size_t m = 0; m < md.size(); ++m) {
                double v = md[m];
                if (!std::isfinite(v)) return false;
                double d = std::sqrt(v);
                double w = magnetWeight[m];
                weighted_sum_d2 += v * w;
                weighted_sum_d  += d * w;
                sum_d += d;
                if (d > max_d) max_d = d;
                if (outResidual) (*outResidual)[i * M + m] = std::sqrt(w) * d;
                total++;
            }
        }
        // ============== ПОЗИЦИОННАЯ КОМПОНЕНТА (БЕЗ ВЕСА) ==============
        // position-loss активна, если cfg.usePositionLoss (в "envelope" mode — выключена).
        double pos_unweighted = cfg.usePositionLoss
                      ? (cfg.lossSquared ? weighted_sum_d2 : weighted_sum_d)
                      : 0.0;
        if (outPosComp) *outPosComp = pos_unweighted;
        outLoss = pos_unweighted;  // position имеет фактический вес 1.0
        // Max-штраф (давит худшую точку). max_d в метрах; max_d² ~ pos (м²).
        if (cfg.lossMaxWeight > 0.0 && cfg.usePositionLoss)
            outLoss += cfg.lossMaxWeight * (cfg.lossSquared ? (max_d*max_d) : max_d);
        // Меж-частичный член: Σ_магнит Σ_{i<j} |closePos[i][m]-closePos[j][m]|².
        // Штрафует расхождение пучка (отделение частиц). Тот же масштаб, что pos.
        if (interOn) {
            double inter = 0.0;
            for (size_t m = 0; m < M; ++m)
                for (size_t i = 0; i < beam.size(); ++i)
                    for (size_t j = i + 1; j < beam.size(); ++j) {
                        double dx = closePos[i][m].x - closePos[j][m].x;
                        double dy = closePos[i][m].y - closePos[j][m].y;
                        double dz = closePos[i][m].z - closePos[j][m].z;
                        inter += dx*dx + dy*dy + dz*dz;
                    }
            outLoss += cfg.interparticleWeight * inter;
        }
        // ОГИБАЮЩАЯ через НАКЛОН ЛИНИИ ТРЕНДА (а не сумму превышений над стартом):
        // a_m = sqrt(Σ_{i<j}|p_i-p_j|²) -- амплитуда взаимного разлёта частиц в точке m.
        // Линейный МНК-фит a_m ≈ α + β·x, где x = m/(M-1) ∈ [0,1]. β -- подъём тренда
        // от начала к концу окна (в метрах). Штраф = max(0, β): огибающая должна НЕ
        // ВОЗРАСТАТЬ (β ≤ 0). Осцилляции усредняются фитом -- ловится чистый рост, а не
        // транзиентные пики. Масштаб (м) не зависит от числа точек M.
        if (decayOn) {
            double sx = 0.0, sy = 0.0, sxy = 0.0, sxx = 0.0;
            size_t Mc = 0;
            for (size_t m = 0; m < M; ++m) {
                double Sm = 0.0;
                for (size_t i = 0; i < beam.size(); ++i)
                    for (size_t j = i + 1; j < beam.size(); ++j) {
                        double dx = closePos[i][m].x - closePos[j][m].x;
                        double dy = closePos[i][m].y - closePos[j][m].y;
                        double dz = closePos[i][m].z - closePos[j][m].z;
                        Sm += dx*dx + dy*dy + dz*dz;
                    }
                double a = std::sqrt(Sm);
                double x = (M > 1) ? (double)m / (double)(M - 1) : 0.0;
                sx += x; sy += a; sxy += x*a; sxx += x*x; ++Mc;
            }
            double slope = 0.0;
            if (Mc >= 2) {
                double denom = (double)Mc * sxx - sx * sx;
                if (std::abs(denom) > 1e-30)
                    slope = ((double)Mc * sxy - sx * sy) / denom;
            }
            double penalty = (slope > 0.0) ? slope : 0.0;   // штраф только за РОСТ
            if (outDecayComp) *outDecayComp = penalty;
            outLoss += cfg.envelopeDecayWeight * penalty;
        }
        // ПАРАЛЛЕЛЬНОСТЬ движения: спред направлений импульса частиц на финише
        // (тот же член, что в оконном loss). parallel = Σ_{i<j}|p̂_i - p̂_j|².
        // Штрафует угловую расходимость; среди прочего давит чит-расщепление
        // (у разделившихся групп импульсы расходятся -> большой штраф).
        if (cfg.parallelismWeight > 0.0) {
            std::vector<FieldVector> pdir(beam.size());
            for (size_t i = 0; i < beam.size(); ++i) {
                const auto& mom = beam[i].mom;
                double pn = std::sqrt(mom.x*mom.x + mom.y*mom.y + mom.z*mom.z);
                if (pn < 1e-30) pn = 1e-30;
                pdir[i] = {mom.x/pn, mom.y/pn, mom.z/pn};
            }
            double parallel = 0.0;
            for (size_t i = 0; i < beam.size(); ++i)
                for (size_t j = i + 1; j < beam.size(); ++j) {
                    double dx = pdir[i].x - pdir[j].x;
                    double dy = pdir[i].y - pdir[j].y;
                    double dz = pdir[i].z - pdir[j].z;
                    parallel += dx*dx + dy*dy + dz*dz;
                }
            if (outParallelComp) *outParallelComp = parallel;
            outLoss += cfg.parallelismWeight * parallel;
        }
        // mean/max-mm считаем всегда (они нужны для диагностики, не loss).
        outMeanMM = (total > 0) ? (sum_d / total) * 1000.0 : 0.0;
        outMaxMM  = max_d * 1000.0;

        // ============== ПЕРИОДИЧЕСКИЙ ЧЛЕН: выход ЦЕНТРОИДА ∥ орбите ==============
        // Замкнутая периодическая орбита: ЦЕНТРОИД пучка (частица 0, стартует на орбите)
        // должен выходить из супериода параллельно касательной орбиты fpDir — тогда SP2,
        // SP3… повторяют SP1. Спутники НЕ форсируем (они трассируют бетатрон; их «матчинг»
        // держат parallelism/envelope). Член добавляется И в loss (видит line_search,
        // устойчивый к джиттеру), И в хвост невязки SVD: w²·(a1²+a2²), где a1,a2 — поперечные
        // (⊥ fpDir) направляющие косинусы импульса центроида.
        if (angleActive) {
            double nn = std::sqrt(fpDir.x*fpDir.x + fpDir.y*fpDir.y + fpDir.z*fpDir.z);
            FieldVector nhat{fpDir.x/nn, fpDir.y/nn, fpDir.z/nn};
            FieldVector e1{-nhat.y, nhat.x, 0.0};            // поперечная в плоскости орбиты
            double e1n = std::sqrt(e1.x*e1.x + e1.y*e1.y + e1.z*e1.z);
            if (e1n < 1e-9) { e1 = {1.0,0.0,0.0}; e1n = 1.0; }   // вырожд.: nhat вертикальна
            e1 = {e1.x/e1n, e1.y/e1n, e1.z/e1n};
            FieldVector e2{ nhat.y*e1.z - nhat.z*e1.y,        // e2 = nhat × e1 (ортонормир.)
                            nhat.z*e1.x - nhat.x*e1.z,
                            nhat.x*e1.y - nhat.y*e1.x };
            const double wA = cfg.svdAngleWeight;
            const auto& p = beam[0].mom;                      // ЦЕНТРОИД
            double pn = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
            if (pn < 1e-30) pn = 1e-30;
            double a1 = (p.x*e1.x + p.y*e1.y + p.z*e1.z) / pn;
            double a2 = (p.x*e2.x + p.y*e2.y + p.z*e2.z) / pn;
            outLoss += wA*wA * (a1*a1 + a2*a2);               // в composite-loss (line_search)
            if (angleRes) {                                   // в невязку SVD
                (*outResidual)[beam.size()*M]     = wA * a1;
                (*outResidual)[beam.size()*M + 1] = wA * a2;
            }
        }

        // ============== УГЛОВАЯ КОМПОНЕНТА LOSS ==============
        // Σ_по_частицам sin²(θ_i), без веса. Возвращается через outDirComp,
        // в outLoss входит с базовым весом cfg.dirWeight.
        double dir_unweighted = 0.0;
        if (cfg.dirWeight > 0.0 || outDirComp) {
            // Считаем если: (a) включён базовым весом, (b) запрошен dynamic-режимом.
            for (size_t i = 0; i < beam.size(); ++i) {
                const auto& p = beam[i].mom;
                double pn = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
                if (pn < 1e-30) continue;
                double dot = (p.x * cfg.dirReference.x +
                              p.y * cfg.dirReference.y +
                              p.z * cfg.dirReference.z) / pn;
                if (dot > 1.0) dot = 1.0;
                if (dot < -1.0) dot = -1.0;
                dir_unweighted += 1.0 - dot * dot;
            }
            if (cfg.dirWeight > 0.0) outLoss += cfg.dirWeight * dir_unweighted;
        }
        if (outDirComp) *outDirComp = dir_unweighted;

        // ============== ENVELOPE-РОСТ (отклонение от ДИЗАЙН-ОРБИТЫ) ==============
        // Старая формулировка использовала P1 как референс, но P1 сама может
        // уплывать -- тогда штраф несправедлив (показывает «спутник ушёл от P1»,
        // хотя на самом деле это P1 ушла от орбиты).
        //
        // Новая: для КАЖДОЙ частицы (включая P1!) считаем расстояние до
        // полилинии дизайн-орбиты:
        //   A_init_i  = distance(beam0[i].pos, designPolyline)
        //   A_final_i = distance(beam[i].pos,  designPolyline)
        //   growth_i  = max(0, A_final_i - A_init_i)
        // Loss_env = envWeight × Σ_i growth_i².
        //
        // Для P1 (стартует на орбите) A_init ≈ 0, и любой drift даёт штраф.
        // Для спутников (старт ±1 мм) штраф только за РОСТ амплитуды.
        double env_unweighted = 0.0;
        if ((cfg.envWeight > 0.0 || outEnvComp) && designPolyline.size() >= 2) {
            for (size_t i = 0; i < beam.size(); ++i) {
                double a_init  = distanceToPolyline(beam0[i].pos, designPolyline);
                double a_final = distanceToPolyline(beam[i].pos,  designPolyline);
                double growth = a_final - a_init;
                if (growth > 0.0) env_unweighted += growth * growth;
            }
            if (cfg.envWeight > 0.0) outLoss += cfg.envWeight * env_unweighted;
        }
        if (outEnvComp) *outEnvComp = env_unweighted;

        return true;
    };

    // Старая runTracking -- используется ТОЛЬКО для финального вывода .bin
    // (в нём всё ещё хранится полная траектория, чтобы visualize.py мог отрисовать).
    auto runTracking = [&](std::vector<std::vector<FieldVector>>& results) -> bool {
        // reset particles
        for (size_t i = 0; i < beam.size(); ++i) beam[i] = beam0[i];
        results.assign(beam.size(), {});
        // FinishPlane хранит callCount per-instance -- по копии на частицу.
        std::vector<FinishPlane> planes(beam.size());
        for (size_t i = 0; i < beam.size(); ++i) {
            planes[i].setup(fp1, fp2, fp3, startRef, maxTurns, fpDir);
        }
        bool ok = true;
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < (int)beam.size(); ++i) {
            int steps = 0;
            results[i].reserve(20000);
            results[i].push_back(beam[i].pos);  // стартовая точка
            while (steps < 50000000) {
                Integrator::borisStep(beam[i], masses[i], charges[i], lattice, dt);
                steps++;
                bool crossed = planes[i].hasCrossed(beam[i].pos);
                if ((steps % trajStride) == 0 || crossed) {
                    results[i].push_back(beam[i].pos);
                }
                if (crossed) break;
            }
        }
        // Проверка валидности: все частицы должны успеть пройти достаточно
        for (size_t i = 0; i < results.size(); ++i)
            if (results[i].size() < 2) { ok = false; break; }
        return ok;
    };

    // ПРИМЕЧАНИЕ: Старая computeLoss с воксельной сеткой удалена.
    // Теперь loss считается ИНКРЕМЕНТАЛЬНО прямо в runTrackingAndLoss()
    // (см. выше). Эффект: 1) не нужна траектория в памяти, 2) пропускается
    // post-tracking-проход по сетке для kdtree-like поиска.

    // ---- 10. Baseline ----
    std::vector<std::vector<FieldVector>> traj;  // используется только для финального .bin
    double curLoss = 0, curMeanMM = 0, curMaxMM = 0;
    // Компоненты loss БЕЗ весов -- нужны для dynamic_loss reweighting (MC).
    // Сохраняются как cur_pos/cur_dir/cur_env, обновляются при accept.
    double cur_pos_c = 0, cur_dir_c = 0, cur_env_c = 0;
    if (!runTrackingAndLoss(curLoss, curMeanMM, curMaxMM,
                             &cur_pos_c, &cur_dir_c, &cur_env_c)) {
        std::cerr << "ERROR: baseline tracking failed\n";
        return 1;
    }
    double bestLoss = curLoss;
    std::vector<double> bestValues(rgroups.size());
    for (size_t i = 0; i < rgroups.size(); ++i) bestValues[i] = rgroups[i].currentValue;
    std::cout << "\nBaseline: loss=" << std::fixed << std::setprecision(6) << curLoss
              << "  mean=" << std::setprecision(2) << curMeanMM << "mm"
              << "  max="  << curMaxMM << "mm\n\n";

    // ---- 11. History accumulator ----
    struct HistEntry {
        int it;
        std::string name;
        double oldVal, newVal;
        double loss, bestLoss;
        double meanMM, maxMM;
        std::string mark;
    };
    std::vector<HistEntry> history;
    history.reserve(cfg.nIter + 16);

    // ---- 12. Optimizer loop ----
    std::mt19937 rng(cfg.seed);
    std::uniform_real_distribution<double> uniMP(-1.0, 1.0);
    int accepted = 0, rejected = 0, kicks = 0, failed = 0, reverts = 0;
    int normalIdx = 0;
    int itersSinceKick = 0;
    int kickTrigger = cfg.kickEvery;
    int kicksWithoutBest = 0;       // подряд кикcов без нового best
    const int REVERT_AFTER_KICKS = cfg.revertAfterKicks;  // 0 = disable

    auto randomStep = [&](RuntimeGroup& rg, double scale) -> double {
        if (rg.cfg.stepAbs > 0)
            return rg.cfg.stepAbs * scale * uniMP(rng);
        double base = std::abs(rg.currentValue);
        if (base < 1e-12) base = 1.0;
        return base * rg.cfg.stepPercent * scale * uniMP(rng);
    };

    StepTimer loopTimer;
    // Сохранение истории в JSON. Пишется (а) после каждых CHECKPOINT_EVERY итераций
    // как чекпойнт (на случай Ctrl+C / краша), (б) обязательно в конце.
    // Для чекпойнтов используем БЕСТ-значения из bestValues (а не текущие),
    // чтобы при прерывании можно было сохранить best в CSV.
    const std::string histPath = ws + "results/" + cfg.historyOut;
    auto writeHistory = [&](bool finalCall) {
        std::ofstream hf(histPath);
        if (!hf.is_open()) {
            std::cerr << "WARN: cannot write " << histPath << "\n";
            return;
        }
        double elapsed = loopTimer.sec();
        hf << "{\n";
        hf << "  \"experiment_id\": " << cfg.expId << ",\n";
        hf << "  \"n_iterations\": " << cfg.nIter << ",\n";
        hf << "  \"completed_iter\": " << (int)history.size() << ",\n";
        hf << "  \"final\": " << (finalCall ? "true" : "false") << ",\n";
        hf << "  \"best_loss\": "    << std::fixed << std::setprecision(10) << bestLoss << ",\n";
        hf << "  \"accepted\": "  << accepted << ",\n";
        hf << "  \"rejected\": "  << rejected << ",\n";
        hf << "  \"kicks\": "     << kicks << ",\n";
        hf << "  \"reverts\": "   << reverts << ",\n";
        hf << "  \"failed\": "    << failed << ",\n";
        hf << "  \"loop_sec\": "  << std::fixed << std::setprecision(3) << elapsed << ",\n";
        hf << "  \"group_names\": [";
        for (size_t k = 0; k < rgroups.size(); ++k)
            hf << (k ? "," : "") << "\"" << rgroups[k].cfg.name << "\"";
        hf << "],\n";
        hf << "  \"best_values\": [";
        hf << std::setprecision(10);
        for (size_t k = 0; k < bestValues.size(); ++k)
            hf << (k ? "," : "") << bestValues[k];
        hf << "],\n";
        // Полное описание групп (для apply_optim_history.py).
        hf << "  \"groups\": [\n";
        for (size_t k = 0; k < rgroups.size(); ++k) {
            const auto& rg = rgroups[k];
            hf << "    {\"name\":\"" << rg.cfg.name << "\""
               << ",\"target\":\"" << rg.cfg.target << "\""
               << ",\"file\":\""   << rg.cfg.file   << "\""
               << ",\"value_col\":\"" << rg.cfg.valueCol << "\""
               << ",\"best\":" << std::setprecision(10) << bestValues[k]
               << ",\"matchers\":[";
            for (size_t mi = 0; mi < rg.cfg.matchers.size(); ++mi) {
                const auto& m = rg.cfg.matchers[mi];
                hf << (mi ? "," : "") << "{\"column\":\"" << m.column << "\",\"values\":[";
                for (size_t vi = 0; vi < m.values.size(); ++vi)
                    hf << (vi ? "," : "") << "\"" << m.values[vi] << "\"";
                hf << "]}";
            }
            hf << "]}";
            if (k + 1 < rgroups.size()) hf << ",";
            hf << "\n";
        }
        hf << "  ],\n";
        hf << "  \"iterations\": [\n";
        for (size_t k = 0; k < history.size(); ++k) {
            const auto& e = history[k];
            hf << "    {\"it\":" << e.it
               << ",\"name\":\"" << e.name << "\""
               << ",\"old\":" << e.oldVal
               << ",\"new\":" << e.newVal
               << ",\"loss\":" << e.loss
               << ",\"best\":" << e.bestLoss
               << ",\"mean_mm\":" << std::setprecision(3) << e.meanMM
               << ",\"max_mm\":" << e.maxMM
               << std::setprecision(10)
               << ",\"mark\":\"" << e.mark << "\"}";
            if (k + 1 < history.size()) hf << ",";
            hf << "\n";
        }
        hf << "  ]\n";
        hf << "}\n";
    };
    const int CHECKPOINT_EVERY = 100;

    // ========================================================================
    // PER-ELEMENT MODE — скользящее окно по магнитам.
    // На каждом положении окна:
    //   - Стартуем 5 частиц НА ОПОРНОЙ ОРБИТЕ перед окном (P1 на ней, спутники
    //     ±5 мм поперёк tangent в локальной системе).
    //   - Симулируем только через активные магниты окна (mask по global_id).
    //   - Финиш-плоскость = перпендикуляр к опорной орбите за выходом окна.
    //   - m итераций hill-climber MC, варьируя ТОЛЬКО группы окна.
    //   - Loss = ||P1.exit - target||² + dirWeight × (1 - dot(P1_tangent, target_tangent))².
    // Окно сдвигается с шагом stride, N_passes полных проходов по решётке.
    // ========================================================================
    if (cfg.method == "windowed" || cfg.method == "lens") {
        // LENS-режим: окно скользит по ЛИНЗАМ, в начало окна инжектируется свежий
        // прямой пучок (вместо prefix-snapshot сжатого пучка от старта конфигурации).
        const bool lensWindowing = (cfg.method == "lens")
                                || (cfg.perElement.windowUnit == "lens");
        if (lensWindowing) {
            if (cfg.perElement.lensDeviceTypes.empty()) {
                cfg.perElement.lensDeviceTypes.insert("QP");
                cfg.perElement.lensDeviceTypes.insert("DQ");
            }
            // Инжекция свежего пучка в каждое окно несовместима с
            // prefix-from-start / accumulate / lensLag / drift-фазой — отключаем.
            cfg.perElement.startFromSpStart = false;
            cfg.perElement.accumulateGroups = false;
            cfg.perElement.lensLagWindows   = 0;
            cfg.perElement.driftExtensionM  = 0.0;
            cfg.perElement.groupDipoleArrays = false;
            std::string lt;
            for (const auto& s : cfg.perElement.lensDeviceTypes)
                lt += (lt.empty() ? "" : ",") + s;
            std::cout << "\n[LENS-РЕЖИМ (окно по линзам, свежий прямой пучок)]"
                      << "  линз_в_окне=" << cfg.perElement.lensWindowSize
                      << "  линзы={" << lt << "}"
                      << "  спред=±" << cfg.perElement.injectSpreadM * 1000.0 << "мм\n"
                      << "  loss = position + parallelism + envelope_decay (локально в окне);"
                      << " первое окно от начала конфигурации, последнее до глобальной плоскости\n";
        } else {
            std::cout << "\n[ОКОННЫЙ РЕЖИМ (windowed)]"
                      << " W=" << cfg.perElement.windowSize
                      << " stride=" << cfg.perElement.windowStride
                      << " m=" << cfg.perElement.itersPerWindow
                      << " passes=" << cfg.perElement.passes << "\n";
        }
        // Явная диагностика режима поиска внутри окна.
        if (cfg.perElement.searchMode == "line_search") {
            std::cout << "  search_mode=LINE_SEARCH (последовательное приближение)"
                      << "  K=" << std::max(3, cfg.perElement.lineSearchPoints)
                      << "  rounds<=" << cfg.perElement.lineSearchMaxRounds
                      << "  диапазоны: |ток|<=" << cfg.perElement.currentRangeAmp << "A"
                      << ", |coeff|<=" << cfg.perElement.coeffRange << "\n";
            if (cfg.perElement.lineSearchLocalCurrent > 0.0 || cfg.perElement.lineSearchLocalCoeff > 0.0)
                std::cout << "  ЛОКАЛЬНЫЙ поиск (в окне): окрестность вокруг текущего  ток ±"
                          << cfg.perElement.lineSearchLocalCurrent << "A, coeff ±"
                          << cfg.perElement.lineSearchLocalCoeff << "  (против railing'а)\n";
        } else {
            std::cout << "  search_mode=MONTE_CARLO (случайные кики — метод Монте-Карло)\n";
        }
        if (cfg.perElement.lossMode == "composite_next") {
            if (!cfg.pathPoints.empty())
                std::cout << "  loss=composite_next, ОПОРНЫЙ ПУТЬ: " << cfg.pathPoints.size()
                          << " точек (откл. от пути, а не от центров)\n";
            else
                std::cout << "  loss=composite_next: PATHPT не заданы -> fallback на центры магнитов\n";
        }

        // 1. Загрузить magnet_endpoints CSV
        if (cfg.magnetEndpointsCsv.empty()) {
            std::cerr << "ERROR: оконный режим (windowed) requires STRING MAGNET_ENDPOINTS_CSV "
                         "(сгенерируйте через prepare_optim_config.py)\n";
            return 1;
        }
        auto endpoints = loadMagnetEndpointsCsv(cfg.magnetEndpointsCsv);
        if (endpoints.empty()) {
            std::cerr << "ERROR: magnet endpoints не загружены ("
                      << cfg.magnetEndpointsCsv << ")\n";
            return 1;
        }
        // Индексируем endpoints по global_id для быстрого lookup
        std::unordered_map<int, const MagnetEndpoint*> endpointsByGid;
        for (const auto& e : endpoints) endpointsByGid[e.global_id] = &e;
        std::cout << "  loaded " << endpoints.size() << " magnet endpoints\n";

        // 2. Построить map global_id -> [group_index_in_rgroups]
        // Группа регистрируется ПО ВСЕМ её gid'ам. Это значит: для lattice
        // с несколькими SP одна per-position группа покрывает несколько gid
        // (например, QF1A в каждом из 4 SP). Она будет «активна» в каждом
        // окне, куда попал хотя бы один её gid. Мутация изменит ВСЕ
        // instance-ы (через applyGroupValue), но симуляция окна видит
        // только активные через mask магниты -- то есть эффект «one knob
        // per type» сохраняется глобально, а локально -- нормальный per_element.
        std::unordered_map<int, std::vector<int>> groupsByGid;
        int n_multi_gid = 0;
        for (int gi = 0; gi < (int)rgroups.size(); ++gi) {
            const auto& rg = rgroups[gi];
            if (rg.elementIndices.empty()) continue;
            std::set<int> gids;
            for (int ei : rg.elementIndices) gids.insert(meta[ei].global_id);
            if (gids.size() > 1) n_multi_gid++;
            for (int gid : gids) groupsByGid[gid].push_back(gi);
        }
        if (n_multi_gid > 0)
            std::cout << "  note: " << n_multi_gid
                      << " групп покрывают несколько gid (multi-SP) -- "
                         "регистрирую в каждом\n";

        // 3. Список активных gids: пересечение endpoints с груповыми
        std::vector<int> activeGids;
        for (const auto& e : endpoints)
            if (groupsByGid.count(e.global_id)) activeGids.push_back(e.global_id);
        std::sort(activeGids.begin(), activeGids.end());
        std::cout << "  активных магнитов: " << activeGids.size() << "\n";

        const int W      = lensWindowing ? std::max(1, cfg.perElement.lensWindowSize)
                                         : cfg.perElement.windowSize;
        const int stride = std::max(1, cfg.perElement.windowStride);
        const int m_iter = cfg.perElement.itersPerWindow;
        const int passes = cfg.perElement.passes;
        if (!lensWindowing && (W < 1 || (int)activeGids.size() < W)) {
            std::cerr << "ERROR: окно W=" << W << " больше числа активных магнитов "
                      << activeGids.size() << "\n";
            return 1;
        }

        // 4. Базовая скорость частицы (для спавна на эталоне) -- из P1 в beam0.
        // Все 5 частиц получают эту же |v| (релятивистская: |v| < c).
        auto momentumToVelocity = [](const FieldVector& p, double m0) -> FieldVector {
            double p2 = p.x*p.x + p.y*p.y + p.z*p.z;
            double gamma = std::sqrt(1.0 + p2 / (m0 * m0 * C_SQ));
            return p * (1.0 / (m0 * gamma));
        };
        if (beam0.empty()) {
            std::cerr << "ERROR: beam0 пуст -- не из чего взять базовую скорость\n";
            return 1;
        }
        FieldVector v0_full = momentumToVelocity(beam0[0].mom, masses[0]);
        const double v0_mag = std::sqrt(v0_full.x*v0_full.x +
                                        v0_full.y*v0_full.y +
                                        v0_full.z*v0_full.z);
        const double satellite_offset_m = 0.005;  // 5 мм
        const int MAX_STEPS_PER_WINDOW  = 200000;
        // Direction-вес для per_element loss (если задан в JSON через
        // direction_loss.weight -- используем его; иначе 0).
        const double dirW = cfg.dirWeight;

        // 5. Lambda: одна симуляция окна -> возвращает loss или infinity при сбое.
        // Использует текущее состояние lattice (с активной маской и параметрами).
        // prefixSnapshot (опц.): если != nullptr и непустой, частицы стартуют
        // ИЗ ЭТОГО состояния (готовый "реальный пучок" на entry окна) вместо
        // spawn'a на эталонной орбите. Это даёт оптимизацию prefix-run.
        auto runWindowSim = [&](const FieldVector& start_pos, double start_yaw,
                                const FieldVector& target_pt, double target_yaw,
                                const std::vector<ParticleState>* prefixSnapshot = nullptr) -> double {
            // Эталонные единичные вектора
            FieldVector t_start  = {std::cos(start_yaw),  std::sin(start_yaw),  0.0};
            FieldVector y_loc    = {-std::sin(start_yaw), std::cos(start_yaw),  0.0};
            FieldVector t_target = {std::cos(target_yaw), std::sin(target_yaw), 0.0};
            FieldVector y_tgt    = {-std::sin(target_yaw),std::cos(target_yaw), 0.0};
            FieldVector z_axis   = {0.0, 0.0, 1.0};

            // Стартовое состояние частиц: либо snapshot (если задан), либо
            // spawn на эталонной орбите.
            const int NP = 5;
            std::vector<ParticleState> p(NP);
            std::vector<double>        pm(NP, masses[0]);
            std::vector<double>        pq(NP, charges[0]);
            if (prefixSnapshot && (int)prefixSnapshot->size() >= NP) {
                for (int i = 0; i < NP; ++i) p[i] = (*prefixSnapshot)[i];
            } else {
                FieldVector v_init = t_start * v0_mag;
                p[0].pos = start_pos;
                p[0].mom = Integrator::velocityToMomentum(v_init, masses[0]);
                p[1].pos = {start_pos.x + satellite_offset_m * y_loc.x,
                            start_pos.y + satellite_offset_m * y_loc.y, start_pos.z};
                p[2].pos = {start_pos.x - satellite_offset_m * y_loc.x,
                            start_pos.y - satellite_offset_m * y_loc.y, start_pos.z};
                p[3].pos = {start_pos.x, start_pos.y, start_pos.z + satellite_offset_m};
                p[4].pos = {start_pos.x, start_pos.y, start_pos.z - satellite_offset_m};
                for (int i = 1; i < NP; ++i)
                    p[i].mom = Integrator::velocityToMomentum(v_init, masses[0]);
            }

            // Финиш-плоскость через 3 точки на ней
            FieldVector fp_p0 = target_pt;
            FieldVector fp_p1 = {target_pt.x + y_tgt.x,
                                 target_pt.y + y_tgt.y,
                                 target_pt.z + y_tgt.z};
            FieldVector fp_p2 = {target_pt.x + z_axis.x,
                                 target_pt.y + z_axis.y,
                                 target_pt.z + z_axis.z};

            std::vector<FinishPlane> planes(NP);
            std::vector<char> reached(NP, 0);
            for (int i = 0; i < NP; ++i) planes[i].setup(fp_p0, fp_p1, fp_p2, p[i].pos, 0);

            // path_penalty: накапливается ТОЛЬКО для P1 (центральной), по
            // расстоянию от текущей позиции до design_polyline на каждом
            // pathSampleStride-м шаге пушера.
            const bool    pathOnly = cfg.perElement.pathOnly;
            const double  pathW    = pathOnly ? 1.0 : cfg.perElement.pathWeight;
            const int     pathStr  = std::max(1, cfg.perElement.pathSampleStride);
            const bool    pathOn   = (pathOnly || pathW > 0.0)
                                     && designPolyline.size() >= 2;
            double path_penalty = 0.0;

            // Симуляция последовательно (5 частиц -- слишком мало для OpenMP overhead'а)
            for (int i = 0; i < NP; ++i) {
                int steps = 0;
                while (steps < MAX_STEPS_PER_WINDOW) {
                    Integrator::borisStep(p[i], pm[i], pq[i], lattice, dt);
                    steps++;
                    // path penalty -- только для P1, разреженной выборкой
                    if (pathOn && i == 0 && (steps % pathStr) == 0) {
                        double d = distanceToPolyline(p[0].pos, designPolyline);
                        path_penalty += d * d;
                    }
                    if (planes[i].hasCrossed(p[i].pos)) { reached[i] = 1; break; }
                }
            }
            // Все 5 должны достичь финиша
            for (int i = 0; i < NP; ++i)
                if (!reached[i]) return std::numeric_limits<double>::infinity();

            // Loss
            // ||P1 - target||² -- основная составляющая
            double dx = p[0].pos.x - target_pt.x;
            double dy = p[0].pos.y - target_pt.y;
            double dz = p[0].pos.z - target_pt.z;
            double pos_loss = dx*dx + dy*dy + dz*dz;
            // спутники: ||Pi - target||² (с меньшим весом, как proxy для envelope)
            double sat_loss = 0.0;
            for (int i = 1; i < NP; ++i) {
                double sx = p[i].pos.x - target_pt.x;
                double sy = p[i].pos.y - target_pt.y;
                double sz = p[i].pos.z - target_pt.z;
                sat_loss += sx*sx + sy*sy + sz*sz;
            }
            // Угол P1 vs target tangent: (1 - dot)² штрафует малые отклонения квадратично
            double pn = std::sqrt(p[0].mom.x*p[0].mom.x +
                                  p[0].mom.y*p[0].mom.y +
                                  p[0].mom.z*p[0].mom.z);
            double cos_th = (pn > 1e-30)
                ? (p[0].mom.x * t_target.x + p[0].mom.y * t_target.y + p[0].mom.z * t_target.z) / pn
                : 1.0;
            double dir_loss = (1.0 - cos_th) * (1.0 - cos_th);

            if (pathOnly) {
                // Чистый режим: только сумма квадратов отклонений P1 от опорной орбиты.
                return path_penalty;
            }
            return pos_loss + 0.25 * sat_loss + dirW * dir_loss + pathW * path_penalty;
        };

        // 5b. Lambda: composite-режим для окна. Использует ту же loss-формулу
        // что и глобальный runTrackingAndLoss, но magnetPos фильтруется по
        // gid окна (extra TARGET'ы с gid=-1 всегда включаются). Position
        // считается инкрементально: на каждом шаге пушера обновляется
        // minDist²[particle][magnet_in_window]. Envelope -- по design_polyline
        // (общая). Direction -- угол P1 относительно target_yaw окна.
        auto runWindowSimComposite = [&](const FieldVector& start_pos, double start_yaw,
                                          const FieldVector& target_pt, double target_yaw,
                                          const std::vector<int>& magIdxs,
                                          const std::vector<ParticleState>* prefixSnapshot = nullptr,
                                          double* outMeanMM = nullptr) -> double {
            const FieldVector t_target = {std::cos(target_yaw), std::sin(target_yaw), 0.0};
            const FieldVector y_loc    = {-std::sin(start_yaw), std::cos(start_yaw),  0.0};
            const FieldVector y_tgt    = {-std::sin(target_yaw),std::cos(target_yaw), 0.0};
            const FieldVector z_axis   = {0.0, 0.0, 1.0};

            const int NP = 5;
            std::vector<ParticleState> p(NP), p0(NP);
            if (prefixSnapshot && (int)prefixSnapshot->size() >= NP) {
                for (int i = 0; i < NP; ++i) { p[i] = (*prefixSnapshot)[i]; p0[i] = p[i]; }
            } else {
                FieldVector v_init = {std::cos(start_yaw)*v0_mag,
                                      std::sin(start_yaw)*v0_mag, 0.0};
                p[0].pos = start_pos;
                p[1].pos = {start_pos.x + satellite_offset_m * y_loc.x,
                            start_pos.y + satellite_offset_m * y_loc.y, start_pos.z};
                p[2].pos = {start_pos.x - satellite_offset_m * y_loc.x,
                            start_pos.y - satellite_offset_m * y_loc.y, start_pos.z};
                p[3].pos = {start_pos.x, start_pos.y, start_pos.z + satellite_offset_m};
                p[4].pos = {start_pos.x, start_pos.y, start_pos.z - satellite_offset_m};
                for (int i = 0; i < NP; ++i) {
                    p[i].mom = Integrator::velocityToMomentum(v_init, masses[0]);
                    p0[i] = p[i];
                }
            }

            // Финиш-плоскость окна
            FieldVector fp_p0 = target_pt;
            FieldVector fp_p1 = {target_pt.x + y_tgt.x, target_pt.y + y_tgt.y, target_pt.z};
            FieldVector fp_p2 = {target_pt.x + z_axis.x, target_pt.y + z_axis.y,
                                 target_pt.z + z_axis.z};
            std::vector<FinishPlane> planes(NP);
            std::vector<char> reached(NP, 0);
            for (int i = 0; i < NP; ++i) planes[i].setup(fp_p0, fp_p1, fp_p2, p[i].pos, 0);

            // Инкрементальный position-loss: minDist²[particle][magnet_local_idx]
            const size_t M = magIdxs.size();
            const double INF = std::numeric_limits<double>::infinity();
            std::vector<std::vector<double>> minD2(NP, std::vector<double>(M, INF));
            // Стартовая позиция: обновим начальные min²
            for (int i = 0; i < NP; ++i) {
                for (size_t m = 0; m < M; ++m) {
                    const auto& mp = magnetPos[magIdxs[m]];
                    double dx = p[i].pos.x - mp.x;
                    double dy = p[i].pos.y - mp.y;
                    double dz = p[i].pos.z - mp.z;
                    minD2[i][m] = dx*dx + dy*dy + dz*dz;
                }
            }

            // Симуляция (5 частиц последовательно)
            for (int i = 0; i < NP; ++i) {
                int steps = 0;
                while (steps < MAX_STEPS_PER_WINDOW) {
                    Integrator::borisStep(p[i], masses[0], charges[0], lattice, dt);
                    steps++;
                    // Обновить minDist² ко всем магнитам окна
                    const FieldVector& pp = p[i].pos;
                    for (size_t m = 0; m < M; ++m) {
                        const auto& mp = magnetPos[magIdxs[m]];
                        double dx = pp.x - mp.x;
                        double dy = pp.y - mp.y;
                        double dz = pp.z - mp.z;
                        double d2 = dx*dx + dy*dy + dz*dz;
                        if (d2 < minD2[i][m]) minD2[i][m] = d2;
                    }
                    if (planes[i].hasCrossed(p[i].pos)) { reached[i] = 1; break; }
                }
            }
            for (int i = 0; i < NP; ++i)
                if (!reached[i]) return std::numeric_limits<double>::infinity();

            // Position-loss + gate-mean (по lattice-магнитам gid>0, без extra targets)
            double pos_loss = 0.0;
            double gate_sum_d = 0.0;
            long   gate_cnt = 0;
            for (int i = 0; i < NP; ++i) {
                for (size_t m = 0; m < M; ++m) {
                    double w = magnetWeight[magIdxs[m]];
                    double v = minD2[i][m];
                    if (!std::isfinite(v)) return std::numeric_limits<double>::infinity();
                    pos_loss += (cfg.lossSquared ? v : std::sqrt(v)) * w;
                    int g = (magIdxs[m] < (int)magnetGid.size()) ? magnetGid[magIdxs[m]] : -1;
                    if (g > 0) { gate_sum_d += std::sqrt(v); gate_cnt++; }
                }
            }
            if (outMeanMM)
                *outMeanMM = (gate_cnt > 0) ? (gate_sum_d / gate_cnt) * 1000.0 : 0.0;

            // Direction-loss: (1 - cos θ)² × cfg.dirWeight, где θ -- угол P1 vs target_tangent окна
            double dir_loss = 0.0;
            if (cfg.dirWeight > 0.0) {
                double pn = std::sqrt(p[0].mom.x*p[0].mom.x + p[0].mom.y*p[0].mom.y + p[0].mom.z*p[0].mom.z);
                double cos_th = (pn > 1e-30)
                  ? (p[0].mom.x * t_target.x + p[0].mom.y * t_target.y + p[0].mom.z * t_target.z) / pn
                  : 1.0;
                if (cos_th > 1.0) cos_th = 1.0; if (cos_th < -1.0) cos_th = -1.0;
                dir_loss = (1.0 - cos_th) * (1.0 - cos_th);
            }

            // ENVELOPE-LOSS НАМЕРЕННО ОТКЛЮЧЕНА в PE composite-режимах.
            // Бетатронные колебания спутников -- глобальное свойство связки
            // магнитов, локальное окно не может на неё повлиять.
            // (Раньше: env_loss = Σ_i max(0, A_final-A_init)²; убрано.)

            return pos_loss + cfg.dirWeight * dir_loss;
        };

        // 5b2. Lambda: full-trajectory симуляция от глобального старта (beam0)
        // до ГЛОБАЛЬНОЙ финиш-плоскости эксперимента (из experiments.csv).
        // Не использует mask (все магниты активны), не использует prefix-cache.
        // Используется для loss_mode == "composite_full" -- оптимизатор видит
        // ВЕСЬ путь пучка и финиш на реальной целевой плоскости (например,
        // 50 м после SP1 для exp 4). target_pt окна игнорируется -- финиш и
        // tangent берутся глобальные.
        auto runFullSimComposite = [&](const FieldVector& target_pt,
                                        double           target_yaw,
                                        const std::vector<int>& magIdxs,
                                        double* outMeanMM = nullptr,
                                        double* outPos    = nullptr,
                                        double* outDir    = nullptr,
                                        double* outEnv    = nullptr,
                                        double* outPosAlt = nullptr,
                                        bool    localFinish = false,
                                        const std::vector<FieldVector>* ctrlPosO    = nullptr,
                                        const std::vector<int>*         ctrlGidO    = nullptr,
                                        const std::vector<double>*      ctrlWeightO = nullptr,
                                        // LENS-режим: если injPos!=nullptr, спавн НЕ из beam0, а свежий
                                        // прямой пучок в injPos (центр), импульсы || касательной injYaw,
                                        // спутники ±injSpread в локальной перпенд. плоскости.
                                        const FieldVector* injPos = nullptr,
                                        double injYaw    = 0.0,
                                        double injSpread = 0.0) -> double {
            // Контрольные точки: по умолчанию центры loss-магнитов (magnetPos), но
            // composite_next передаёт ПЛОТНЫЙ опорный путь (pathPoints) -- тогда loss
            // считает отклонение от точек пути, а не от центров.
            const std::vector<FieldVector>& CP = ctrlPosO    ? *ctrlPosO    : magnetPos;
            const std::vector<int>&         CG = ctrlGidO    ? *ctrlGidO    : magnetGid;
            const std::vector<double>&      CW = ctrlWeightO ? *ctrlWeightO : magnetWeight;
            // Альтернативный референс (диагностика iron vs .lte) -- только для
            // magnetPos; при override (путь) альт = тот же путь.
            const std::vector<FieldVector>& CPalt = ctrlPosO ? *ctrlPosO : magnetPosAlt;
            // Direction-target: для глобальной плоскости -- cfg.dirReference (yaw=9°
            // из JSON); для ЛОКАЛЬНОЙ плоскости (composite_next) -- касательная орбиты
            // у следующего магнита (target_yaw), т.е. нить нанизывается вдоль орбиты.
            const FieldVector t_target = localFinish
                ? FieldVector{std::cos(target_yaw), std::sin(target_yaw), 0.0}
                : cfg.dirReference;
            // Диагностический режим: считаем pos по АЛЬТЕРНАТИВНОМУ референсу
            // (magnetPosAlt) + env. Включается только при ненулевых указателях,
            // чтобы не замедлять hot-path оптимизации.
            const bool diag = (outPosAlt != nullptr) || (outEnv != nullptr);

            // Спавн: из beam0 (глобальный старт) ИЛИ свежий прямой пучок в injPos
            // (LENS-режим): центр на орбите, импульсы || касательной injYaw,
            // (NP-1) спутников РАВНОМЕРНО по ОКРУЖНОСТИ радиуса injSpread в локальной
            // перпенд. плоскости (базис y_loc, z). Работает для любого NP.
            const int NP = (int)beam0.size();
            if (NP < 1) return std::numeric_limits<double>::infinity();
            std::vector<ParticleState> p, p0;
            if (injPos) {
                p.assign(NP, ParticleState{});
                const double TWO_PI = 6.283185307179586;
                const FieldVector ylw = {-std::sin(injYaw), std::cos(injYaw), 0.0}; // поперечная в плоскости
                const FieldVector zax = {0.0, 0.0, 1.0};                            // вертикаль
                FieldVector v_init = {std::cos(injYaw) * v0_mag,
                                      std::sin(injYaw) * v0_mag, 0.0};
                FieldVector mom0 = Integrator::velocityToMomentum(v_init, masses[0]);
                p[0].pos = *injPos; p[0].mom = mom0;          // центр
                const int nsat = NP - 1;
                for (int i = 1; i < NP; ++i) {                // спутники по окружности
                    double th = (nsat > 0) ? TWO_PI * (i - 1) / nsat : 0.0;
                    double cy = injSpread * std::cos(th);     // вдоль y_loc
                    double cz = injSpread * std::sin(th);     // вдоль z
                    p[i].pos = {injPos->x + cy*ylw.x + cz*zax.x,
                                injPos->y + cy*ylw.y + cz*zax.y,
                                injPos->z + cy*ylw.z + cz*zax.z};
                    p[i].mom = mom0;
                }
                p0 = p;
            } else {
                p = beam0; p0 = beam0;
            }

            // Финиш-плоскость: ГЛОБАЛЬНАЯ (fp1/fp2/fp3 из experiments.csv, 50м за SP1)
            // ИЛИ ЛОКАЛЬНАЯ (composite_next): перпендикуляр к орбите у target_pt
            // (центр следующего магнита), ориентация по target_yaw. Тогда пучок
            // нанизывается окно-за-окном на следующий магнит, а не летит сразу 50м.
            std::vector<FinishPlane> planes(NP);
            std::vector<char> reached(NP, 0);
            FieldVector lf_a, lf_b, lf_c;
            if (localFinish) {
                const FieldVector y_tgt = {-std::sin(target_yaw), std::cos(target_yaw), 0.0};
                lf_a = target_pt;
                lf_b = {target_pt.x + y_tgt.x, target_pt.y + y_tgt.y, target_pt.z};
                lf_c = {target_pt.x, target_pt.y, target_pt.z + 1.0};
            }
            // Референс плоскости: для свежего пучка (LENS) -- точка инжекции, иначе
            // глобальный startRef (нужен для корректной стороны "до плоскости").
            const FieldVector sref = injPos ? *injPos : startRef;
            for (int i = 0; i < NP; ++i) {
                if (localFinish) planes[i].setup(lf_a, lf_b, lf_c, sref, maxTurns);
                else             planes[i].setup(fp1, fp2, fp3, sref, maxTurns);
            }

            // Инкрементальный min-dist по подмножеству magIdxs (магниты решётки
            // с gid<=gid_last + extra-цели). Loss суммируется по всему подмножеству;
            // gate-mean -- ЛОКАЛЬНЫЙ: только по loss-магнитам решётки (gid>0), т.е.
            // по уже ПРОЙДЕННОЙ части траектории до конца текущего окна, БЕЗ 50м-целей.
            const size_t M = magIdxs.size();
            const double INF = std::numeric_limits<double>::infinity();
            std::vector<std::vector<double>> minD2(NP, std::vector<double>(M, INF));
            // Альтернативный референс (диагностика) -- только в diag-режиме.
            std::vector<std::vector<double>> minD2alt;
            if (diag) minD2alt.assign(NP, std::vector<double>(M, INF));
            // Позиции ближайшего подхода (для меж-частичного члена) -- только если включён.
            const bool interOn = (cfg.interparticleWeight > 0.0);
            const bool decayOn = (cfg.envelopeDecayWeight > 0.0);
            const bool needClose = interOn || decayOn;
            std::vector<std::vector<FieldVector>> closePos;
            if (needClose) closePos.assign(NP, std::vector<FieldVector>(M));
            for (int i = 0; i < NP; ++i) {
                for (size_t m = 0; m < M; ++m) {
                    const auto& mp = CP[magIdxs[m]];
                    double dx = p[i].pos.x - mp.x, dy = p[i].pos.y - mp.y, dz = p[i].pos.z - mp.z;
                    minD2[i][m] = dx*dx + dy*dy + dz*dz;
                    if (needClose) closePos[i][m] = p[i].pos;
                    if (diag) {
                        const auto& ma = CPalt[magIdxs[m]];
                        double ax = p[i].pos.x - ma.x, ay = p[i].pos.y - ma.y, az = p[i].pos.z - ma.z;
                        minD2alt[i][m] = ax*ax + ay*ay + az*az;
                    }
                }
            }

            // Симуляция (mask должен быть OFF на момент вызова)
            const int MAX_FULL_STEPS = MAX_STEPS_PER_WINDOW * 20;
            for (int i = 0; i < NP; ++i) {
                int steps = 0;
                while (steps < MAX_FULL_STEPS) {
                    Integrator::borisStep(p[i], masses[i], charges[i], lattice, dt);
                    steps++;
                    const FieldVector& pp = p[i].pos;
                    for (size_t m = 0; m < M; ++m) {
                        const auto& mp = CP[magIdxs[m]];
                        double dx = pp.x - mp.x, dy = pp.y - mp.y, dz = pp.z - mp.z;
                        double d2 = dx*dx + dy*dy + dz*dz;
                        if (d2 < minD2[i][m]) { minD2[i][m] = d2; if (needClose) closePos[i][m] = pp; }
                        if (diag) {
                            const auto& ma = CPalt[magIdxs[m]];
                            double ax = pp.x - ma.x, ay = pp.y - ma.y, az = pp.z - ma.z;
                            double a2 = ax*ax + ay*ay + az*az;
                            if (a2 < minD2alt[i][m]) minD2alt[i][m] = a2;
                        }
                    }
                    if (planes[i].hasCrossed(p[i].pos)) { reached[i] = 1; break; }
                }
            }
            for (int i = 0; i < NP; ++i)
                if (!reached[i]) return std::numeric_limits<double>::infinity();

            // Position-loss (всё подмножество) + ЛОКАЛЬНЫЙ gate-mean (gid>0)
            // + (diag) pos по альтернативному референсу + max-отклонение².
            double pos_loss = 0.0;
            double pos_alt  = 0.0;
            double maxD2    = 0.0;   // макс. отклонение² по (частица,магнит)
            double gate_sum_d = 0.0; long gate_cnt = 0;
            for (size_t m = 0; m < M; ++m) {
                int absm = magIdxs[m];
                double w = CW[absm];
                bool isLattice = (absm < (int)CG.size()) && (CG[absm] > 0);
                for (int i = 0; i < NP; ++i) {
                    double v = minD2[i][m];
                    if (!std::isfinite(v)) return std::numeric_limits<double>::infinity();
                    pos_loss += (cfg.lossSquared ? v : std::sqrt(v)) * w;
                    if (v > maxD2) maxD2 = v;
                    if (isLattice) { gate_sum_d += std::sqrt(v); gate_cnt++; }
                    if (diag) {
                        double va = minD2alt[i][m];
                        pos_alt += (cfg.lossSquared ? va : std::sqrt(va)) * w;
                    }
                }
            }
            // Max-штраф: давит худшую точку (пик отклонения), а не только сумму.
            double max_pen = (cfg.lossMaxWeight > 0.0)
                ? cfg.lossMaxWeight * (cfg.lossSquared ? maxD2 : std::sqrt(maxD2))
                : 0.0;
            // gate-mean: среднее отклонение пучка от опорной орбиты по уже
            // пройденной части траектории (loss-магниты решётки до gid_last),
            // мм. Без extra 50м-целей -- их недостижимость не должна держать окно.
            if (outMeanMM)
                *outMeanMM = (gate_cnt > 0) ? (gate_sum_d / gate_cnt) * 1000.0 : 0.0;

            // Direction-loss: угол P1 vs target_tangent окна
            double dir_loss = 0.0;
            if (cfg.dirWeight > 0.0) {
                double pn = std::sqrt(p[0].mom.x*p[0].mom.x +
                                      p[0].mom.y*p[0].mom.y +
                                      p[0].mom.z*p[0].mom.z);
                double cos_th = (pn > 1e-30)
                  ? (p[0].mom.x * t_target.x + p[0].mom.y * t_target.y + p[0].mom.z * t_target.z) / pn
                  : 1.0;
                if (cos_th > 1.0) cos_th = 1.0;
                if (cos_th < -1.0) cos_th = -1.0;
                dir_loss = (1.0 - cos_th) * (1.0 - cos_th);
            }

            // ENVELOPE-РОСТ: в оптимизируемый loss НЕ входит (выключена в PE
            // composite -- бетатрон глобален, окно на него не влияет), но в
            // diag-режиме считаем для вывода (та же формула, что в global loss:
            // Σ_i max(0, A_final - A_init)², A = dist до design-полилинии).
            double env_unweighted = 0.0;
            if (diag && designPolyline.size() >= 2) {
                for (int i = 0; i < NP; ++i) {
                    double a_init  = distanceToPolyline(beam0[i].pos, designPolyline);
                    double a_final = distanceToPolyline(p[i].pos,     designPolyline);
                    double growth = a_final - a_init;
                    if (growth > 0.0) env_unweighted += growth * growth;
                }
            }

            // Меж-частичный член: Σ_магнит Σ_{i<j} |closePos[i][m]-closePos[j][m]|².
            // closePos[i][m] = позиция частицы i в точке её ближайшего подхода к
            // контрольной точке m. Тот же масштаб, что pos (per-control-point),
            // штрафует РАСХОЖДЕНИЕ пучка (отделение частиц друг от друга).
            double inter = 0.0;
            if (interOn) {
                for (size_t m = 0; m < M; ++m)
                    for (int i = 0; i < NP; ++i)
                        for (int j = i + 1; j < NP; ++j) {
                            double dx = closePos[i][m].x - closePos[j][m].x;
                            double dy = closePos[i][m].y - closePos[j][m].y;
                            double dz = closePos[i][m].z - closePos[j][m].z;
                            inter += dx*dx + dy*dy + dz*dz;
                        }
            }

            // ============== ОГИБАЮЩАЯ через НАКЛОН ЛИНИИ ТРЕНДА ==============
            // a_m = sqrt(Σ_{i<j}|closePos[i][m]-closePos[j][m]|²) -- амплитуда взаимного
            // разлёта частиц в точке m. Линейный МНК-фит a_m ≈ α + β·x, x=m/(M-1)∈[0,1].
            // decay = max(0, β) -- штраф за РОСТ тренда амплитуды по окну (огибающая
            // должна не возрастать). Осцилляции усредняются фитом; масштаб (м) не
            // зависит от числа точек M.
            double decay = 0.0;
            if (decayOn) {
                double sx = 0.0, sy = 0.0, sxy = 0.0, sxx = 0.0;
                size_t Mc = 0;
                for (size_t m = 0; m < M; ++m) {
                    double Sm = 0.0;
                    for (int i = 0; i < NP; ++i)
                        for (int j = i + 1; j < NP; ++j) {
                            double dx = closePos[i][m].x - closePos[j][m].x;
                            double dy = closePos[i][m].y - closePos[j][m].y;
                            double dz = closePos[i][m].z - closePos[j][m].z;
                            Sm += dx*dx + dy*dy + dz*dz;
                        }
                    double a = std::sqrt(Sm);
                    double x = (M > 1) ? (double)m / (double)(M - 1) : 0.0;
                    sx += x; sy += a; sxy += x*a; sxx += x*x; ++Mc;
                }
                if (Mc >= 2) {
                    double denom = (double)Mc * sxx - sx * sx;
                    double slope = (std::abs(denom) > 1e-30)
                                 ? ((double)Mc * sxy - sx * sy) / denom : 0.0;
                    if (slope > 0.0) decay = slope;
                }
            }

            // ============== ПАРАЛЛЕЛЬНОСТЬ ДВИЖЕНИЯ частиц ==============
            // Спред НАПРАВЛЕНИЙ импульса частиц на финише окна:
            //   parallel = Σ_{i<j} |p̂_i - p̂_j|²,  p̂ = mom/|mom|.
            // Штрафует УГЛОВУЮ расходимость пучка -- частицы должны двигаться
            // параллельно (|p̂_i-p̂_j|² ≈ θ_ij² для малых углов). Дополняет
            // позиционный/огибающий члены: те держат КООРДИНАТЫ вместе, этот -- УГЛЫ.
            double parallel = 0.0;
            if (cfg.parallelismWeight > 0.0) {
                std::vector<FieldVector> pdir(NP);
                for (int i = 0; i < NP; ++i) {
                    const auto& mom = p[i].mom;
                    double pn = std::sqrt(mom.x*mom.x + mom.y*mom.y + mom.z*mom.z);
                    if (pn < 1e-30) pn = 1e-30;
                    pdir[i] = {mom.x/pn, mom.y/pn, mom.z/pn};
                }
                for (int i = 0; i < NP; ++i)
                    for (int j = i + 1; j < NP; ++j) {
                        double dx = pdir[i].x - pdir[j].x;
                        double dy = pdir[i].y - pdir[j].y;
                        double dz = pdir[i].z - pdir[j].z;
                        parallel += dx*dx + dy*dy + dz*dz;
                    }
            }

            // Проброс компонент (диагностика по окну).
            if (outPos)    *outPos    = pos_loss;
            if (outDir)    *outDir    = dir_loss;
            if (outEnv)    *outEnv    = env_unweighted;
            if (outPosAlt) *outPosAlt = pos_alt;

            return pos_loss + cfg.dirWeight * dir_loss + max_pen
                 + cfg.interparticleWeight * inter + cfg.envelopeDecayWeight * decay
                 + cfg.parallelismWeight * parallel;
        };

        // 5c. Lambda: prefix-симуляция от глобального старта до entry окна.
        // Используется при start_from_sp_start=true чтобы получить РЕАЛЬНОЕ
        // состояние пучка на entry окна (с учётом всех мутаций предыдущих
        // окон), сохранить snapshot и переиспользовать его в m итерациях окна.
        //
        // Стартует частицы из beam0 (= particle_groups.csv, обычно (0,0,0)),
        // снимает mask (все магниты активны), запускает каждую частицу до
        // финиш-плоскости на entry окна. Возвращает true если все частицы
        // достигли плоскости. snapshot заполняется состоянием на момент crossing.
        auto runPrefixSim = [&](const FieldVector& entry_pos, double entry_yaw,
                                 std::vector<ParticleState>& snapshot) -> bool {
            lattice.clearActiveMask();  // все магниты активны для prefix
            const int NP = (int)beam0.size();
            if (NP == 0) return false;
            snapshot = beam0;
            // Финиш-плоскость на entry окна (перпендикуляр к tangent)
            FieldVector y_ent = {-std::sin(entry_yaw), std::cos(entry_yaw), 0.0};
            FieldVector z_ax  = {0.0, 0.0, 1.0};
            FieldVector fp0 = entry_pos;
            FieldVector fp1 = {entry_pos.x + y_ent.x, entry_pos.y + y_ent.y, entry_pos.z};
            FieldVector fp2 = {entry_pos.x + z_ax.x,  entry_pos.y + z_ax.y,  entry_pos.z + z_ax.z};
            std::vector<FinishPlane> planes(NP);
            for (int i = 0; i < NP; ++i)
                planes[i].setup(fp0, fp1, fp2, snapshot[i].pos, 0);
            // Лимит шагов: prefix может быть длинным (до 100+ магнитов), даём с запасом.
            const int MAX_PREFIX_STEPS = MAX_STEPS_PER_WINDOW * 20;
            for (int i = 0; i < NP; ++i) {
                int steps = 0;
                bool ok = false;
                while (steps < MAX_PREFIX_STEPS) {
                    Integrator::borisStep(snapshot[i], masses[i], charges[i], lattice, dt);
                    steps++;
                    if (planes[i].hasCrossed(snapshot[i].pos)) { ok = true; break; }
                }
                if (!ok) return false;
            }
            return true;
        };

        // 6. Главный цикл: passes × windows × m итераций MC.
        // После каждого полного pass'а оцениваем ГЛОБАЛЬНЫЙ composite loss
        // (через runTrackingAndLoss -- ту же функцию, что использует обычный MC),
        // чтобы видеть, как локальные per-window улучшения транслируются в
        // глобальную целевую функцию.
        long total_iter = 0;
        StepTimer sw;

        // Глобальная оценка в per_element-режиме: считаем ТОЛЬКО position-loss
        // (= orbit). Direction и envelope бессмысленны, так как локальная
        // настройка магнита не может обеспечить фокусировку -- она зависит
        // от связки магнитов вдоль всей цепочки.
        // Сохраняем оригинальные cfg-флаги, временно ставим orbit-режим,
        // восстанавливаем после каждой оценки.
        auto runGlobalOrbitLoss = [&](double& gL, double& gMean, double& gMax) -> bool {
            bool   sav_usePos = cfg.usePositionLoss;
            double sav_dirW   = cfg.dirWeight;
            double sav_envW   = cfg.envWeight;
            cfg.usePositionLoss = true;
            cfg.dirWeight       = 0.0;
            cfg.envWeight       = 0.0;
            lattice.clearActiveMask();
            bool ok = runTrackingAndLoss(gL, gMean, gMax);
            cfg.usePositionLoss = sav_usePos;
            cfg.dirWeight       = sav_dirW;
            cfg.envWeight       = sav_envW;
            return ok;
        };

        // (Baseline теперь выводится как "[global before pass 1/N]" в начале первого pass'а.)

        // ---- Опционально: map gid → sp_number и sp_number → first gid_in_sp,
        // используется если perElement.startFromSpStart == true ----
        std::unordered_map<int, int> gidToSpNum;
        std::unordered_map<int, int> spNumToFirstGid;
        if (cfg.perElement.startFromSpStart) {
            int g_idx = layoutCsv.col("global_id");
            int s_idx = layoutCsv.col("sp_number");
            if (g_idx < 0 || s_idx < 0) {
                std::cerr << "ERROR: layout has no global_id/sp_number column\n";
                return 1;
            }
            for (const auto& r : layoutCsv.rows) {
                if ((int)r.size() <= std::max(g_idx, s_idx)) continue;
                int g, s;
                try { g = std::stoi(r[g_idx]); s = std::stoi(r[s_idx]); }
                catch (...) { continue; }
                gidToSpNum[g] = s;
                auto it = spNumToFirstGid.find(s);
                if (it == spNumToFirstGid.end() || g < it->second)
                    spNumToFirstGid[s] = g;
            }
            std::cout << "  start_from_sp_start=ON: " << spNumToFirstGid.size()
                      << " superperiods mapped\n";
        }

        // ---- group_dipole_arrays: предвычисление дипольных линеек ----
        // gidDev: device_type по gid; dipArrayEnd[i]: индекс activeGids конца
        // дипольной линейки, к которой принадлежит activeGids[i] (или i, если не
        // диполь / опция выключена). Диполь = device_type начинается с "DP".
        std::unordered_map<int,std::string> gidDev;
        for (const auto& e : meta)
            if (e.global_id > 0 && !gidDev.count(e.global_id))
                gidDev[e.global_id] = e.device_type_id;
        auto isDipGid = [&](int gid) -> bool {
            auto it = gidDev.find(gid);
            return it != gidDev.end() && it->second.rfind("DP", 0) == 0;
        };
        // === LENS-режим: список физических линз (юнитов) для окон ===
        // Линза = magnet device_type из lensDeviceTypes (префикс-матч). Co-located
        // обмотки (DQ_I1/I2/I3 в одной точке) схлопываются в ОДИН юнит, чтобы окно
        // из N линз не вырождалось в одну физическую DQ. repGid = первая обмотка.
        auto isLensGid = [&](int gid) -> bool {
            auto it = gidDev.find(gid);
            if (it == gidDev.end()) return false;
            for (const auto& pref : cfg.perElement.lensDeviceTypes)
                if (it->second.rfind(pref, 0) == 0) return true;
            return false;
        };
        struct LensUnit { int repGid; int gidLo; int gidHi; };
        std::vector<LensUnit> lensUnits;
        if (lensWindowing) {
            for (size_t i = 0; i < endpoints.size(); ) {
                if (!isLensGid(endpoints[i].global_id)) { ++i; continue; }
                const auto& e = endpoints[i];
                size_t j = i;          // схлопнуть co-located обмотки (тот же центр)
                while (j + 1 < endpoints.size()
                       && isLensGid(endpoints[j+1].global_id)
                       && std::abs(endpoints[j+1].x_in  - e.x_in ) < 1e-6
                       && std::abs(endpoints[j+1].y_in  - e.y_in ) < 1e-6
                       && std::abs(endpoints[j+1].x_out - e.x_out) < 1e-6
                       && std::abs(endpoints[j+1].y_out - e.y_out) < 1e-6)
                    ++j;
                lensUnits.push_back({e.global_id, e.global_id, endpoints[j].global_id});
                i = j + 1;
            }
            std::cout << "  линз (юнитов) для окон: " << lensUnits.size() << "\n";
            if ((int)lensUnits.size() < W) {
                std::cerr << "ERROR: линз (" << lensUnits.size()
                          << ") меньше окна W=" << W << "\n";
                return 1;
            }
        }
        const int NA = lensWindowing ? (int)lensUnits.size()
                                     : (int)activeGids.size();
        std::vector<int> dipArrayEnd(lensWindowing ? 0 : NA);
        if (!lensWindowing)
        for (int i = 0; i < NA; ++i) {
            if (!cfg.perElement.groupDipoleArrays || !isDipGid(activeGids[i])) {
                dipArrayEnd[i] = i; continue;
            }
            // Линейка = подряд идущие dipole-slice'ы с ПОСЛЕДОВАТЕЛЬНЫМИ gid
            // (слайсы одного физического диполя, напр. DL1A gid5..9). Это не даёт
            // слить разные диполи, которые в activeGids оказались соседними из-за
            // исключённых промежуточных магнитов (напр. DL1A и DL2B в фазе поворотных).
            int e = i;
            while (e + 1 < NA && isDipGid(activeGids[e + 1])
                   && activeGids[e + 1] == activeGids[e] + 1) ++e;
            dipArrayEnd[i] = e;
        }
        // Конец окна с учётом захвата дипольной линейки целиком.
        auto winEnd = [&](int ws) -> int {
            int we = ws + W - 1;
            if (we >= NA) we = NA - 1;
            if (cfg.perElement.groupDipoleArrays)
                for (int j = ws; j <= we; ++j)
                    if (dipArrayEnd[j] > we) we = dipArrayEnd[j];
            return we;
        };
        // Следующий ws (пропуск нутра захваченной линейки).
        auto nextWs = [&](int ws, int we) -> int {
            int nx = ws + stride;
            if (cfg.perElement.groupDipoleArrays && we > ws + W - 1)
                nx = std::max(nx, we - W + 2);
            return nx;
        };

        // Общее число окон в одном проходе (для нумерации/бюджета).
        int total_windows = 0;
        for (int ws = 0; ws + W <= NA; ) {
            int we = winEnd(ws); ++total_windows; ws = nextWs(ws, we);
        }

        // G = ВСЕГО групп, которые будут активированы (каждая входит "новой"
        // ровно один раз). = все группы по activeGids (windowing их не меняет).
        int G_total = 0;
        {
            std::set<int> allSeen;
            if (lensWindowing) {
                for (const auto& lu : lensUnits)
                    for (int gid = lu.gidLo; gid <= lu.gidHi; ++gid) {
                        auto it = groupsByGid.find(gid);
                        if (it != groupsByGid.end())
                            for (int gi : it->second) allSeen.insert(gi);
                    }
            } else {
                for (int i = 0; i < NA; ++i) {
                    auto it = groupsByGid.find(activeGids[i]);
                    if (it != groupsByGid.end())
                        for (int gi : it->second) allSeen.insert(gi);
                }
            }
            G_total = (int)allSeen.size();
        }
        // ЕДИНЫЙ БЮДЖЕТ: total_iters делится РАВНОМЕРНО по окнам (одна доля на
        // окно), внутри окна line-search распределяет её по группам через
        // budgetRounds() -- тот же алгоритм, что в глобальном режиме и в дрейфе.
        // (Прежний itersPerNew -- "фикс. итер. на новый магнит" -- упразднён;
        // поле осталось лишь как опциональная фаза для random-режима.)
        const int budgetPerWindow =
            (cfg.perElement.totalIters > 0)
              ? std::max(1, cfg.perElement.totalIters / std::max(1, total_windows))
              : 0;
        if (budgetPerWindow > 0)
            std::cout << "  бюджет: total_iters=" << cfg.perElement.totalIters
                      << " РАВНОМЕРНО по " << total_windows << " окнам = "
                      << budgetPerWindow << " симуляций/окно (групп всего G="
                      << G_total << ")\n";

        // Запомним BEFORE-значения для резюме в конце.
        double before_L = 0, before_pos = 0, before_dir = 0, before_env = 0,
               before_mean = 0, before_max = 0;
        bool   before_ok = false;
        for (int pass = 1; pass <= passes; ++pass) {
            // Накопленный набор групп для accumulateGroups режима (сбрасывается на каждом pass'е)
            std::set<int> seenGroups;
            // win_idx, на котором группа ВПЕРВЫЕ появилась в окне (для lagged-lens).
            std::unordered_map<int,int> groupEnterWin;
            // ---- Глобальный composite loss ПЕРЕД проходом ----
            {
                lattice.clearActiveMask();
                double pgL, pgMean, pgMax;
                double pos_c, dir_c, env_c;
                if (runTrackingAndLoss(pgL, pgMean, pgMax, &pos_c, &dir_c, &env_c)) {
                    std::cout << "\n[global before pass " << pass << "/" << passes
                              << "]  composite=" << std::scientific << std::setprecision(3)
                              << pgL
                              << "  pos=" << pos_c
                              << "  dir=" << dir_c
                              << "  env=" << env_c
                              << "  mean=" << std::fixed << std::setprecision(2)
                              << pgMean << "mm  max=" << pgMax << "mm  t="
                              << std::setprecision(1) << totalTimer.sec() << "s\n"
                              << std::defaultfloat;
                    // Запоминаем самое первое (= BEFORE) значение
                    if (pass == 1 && !before_ok) {
                        before_L = pgL; before_pos = pos_c; before_dir = dir_c;
                        before_env = env_c; before_mean = pgMean; before_max = pgMax;
                        before_ok = true;
                    }
                } else {
                    std::cout << "\n[global before pass " << pass << "/" << passes
                              << "]  -- симуляция не дошла до финиша\n";
                }
            }

            int win_idx = 0;
            for (int ws = 0; ws + W <= NA; ) {
                ++win_idx;
                const int we = winEnd(ws);   // конец окна (с захватом дип. линейки)
                // gid_first/gid_last + точка старта/финиша окна.
                int gid_first, gid_last, gid_sim_start, gid_mask_hi;
                FieldVector start_pos, target_pt;
                double start_yaw, target_yaw;
                bool lensFirst = false, lensLast = false, lensLocalFinish = true;
                if (lensWindowing) {
                    // === LENS-режим: окно по линзам, свежий прямой пучок ===
                    lensFirst = (ws == 0);
                    lensLast  = (we == NA - 1);
                    // Первое (от начала конфигурации) и последнее (до глобальной
                    // плоскости) окна содержат на 1 линзу МЕНЬШЕ заявленного размера:
                    // у них есть «граничный» элемент (старт/финиш конфигурации),
                    // занимающий место одной линзы.
                    int lws = ws, lwe = we;
                    if (lensFirst) lwe = std::max(0, std::min(W - 2, NA - 1));
                    if (lensLast)  lws = std::min(NA - 1, std::max(0, (NA - 1) - (W - 2)));
                    if (lws > lwe) lws = lwe;            // защита для крошечных NA
                    const LensUnit& luF = lensUnits[lws];
                    const LensUnit& luL = lensUnits[lwe];
                    gid_first = luF.gidLo;
                    gid_last  = luL.gidHi;          // для групп/скоринга по линзам окна
                    const auto* efp = endpointsByGid[luF.repGid];
                    const auto* elp = endpointsByGid[luL.repGid];
                    if (lensFirst) {
                        // первое окно стартует из НАЧАЛА конфигурации (где реально
                        // стартуют частицы) -> участок входа тоже настраивается.
                        start_pos = beam0.empty() ? FieldVector{0.0,0.0,0.0} : beam0[0].pos;
                        start_yaw = beam0.empty() ? 0.0
                                  : std::atan2(beam0[0].mom.y, beam0[0].mom.x);
                        gid_sim_start = 1;           // от самого начала lattice
                    } else {
                        // свежий пучок инжектируется в СЕРЕДИНУ первой линзы окна.
                        start_pos = {0.5*(efp->x_in+efp->x_out),
                                     0.5*(efp->y_in+efp->y_out), 0.0};
                        start_yaw = std::atan2(efp->y_out - efp->y_in,
                                               efp->x_out - efp->x_in);
                        if (efp->x_in == efp->x_out && efp->y_in == efp->y_out)
                            start_yaw = efp->yaw_in;  // вырожденный сегмент
                        gid_sim_start = luF.gidLo;
                    }
                    if (lensLast) {
                        // последнее окно финиширует на ГЛОБАЛЬНОЙ плоскости конфигурации.
                        target_pt = {0.0, 0.0, 0.0};   // игнор при глобальном финише
                        target_yaw = 0.0;
                        lensLocalFinish = false;
                        gid_mask_hi = std::numeric_limits<int>::max();  // маска до конца
                    } else {
                        // финиш = СЕРЕДИНА последней линзы окна (локальная плоскость).
                        target_pt = {0.5*(elp->x_in+elp->x_out),
                                     0.5*(elp->y_in+elp->y_out), 0.0};
                        target_yaw = std::atan2(elp->y_out - elp->y_in,
                                                elp->x_out - elp->x_in);
                        if (elp->x_in == elp->x_out && elp->y_in == elp->y_out)
                            target_yaw = elp->yaw_in;
                        lensLocalFinish = true;
                        gid_mask_hi = gid_last;
                    }
                } else {
                    // === MAGNET-режим (как раньше): окно по всем активным магнитам ===
                    gid_first = activeGids[ws];
                    gid_last  = activeGids[we];
                    gid_mask_hi = gid_last;
                    const auto* ep_first = endpointsByGid[gid_first];
                    const auto* ep_last  = endpointsByGid[gid_last];
                    // Точка старта: entry первого магнита окна (изолированный старт)
                    // либо entry первого магнита СУПЕРИОДА (накопленный пучок).
                    gid_sim_start = gid_first;
                    start_pos = {ep_first->x_in, ep_first->y_in, 0.0};
                    start_yaw = ep_first->yaw_in;
                    if (cfg.perElement.startFromSpStart) {
                        auto it_sp = gidToSpNum.find(gid_first);
                        if (it_sp != gidToSpNum.end()) {
                            auto it_fg = spNumToFirstGid.find(it_sp->second);
                            if (it_fg != spNumToFirstGid.end()) {
                                int g_sp_start = it_fg->second;
                                auto it_ep = endpointsByGid.find(g_sp_start);
                                if (it_ep != endpointsByGid.end()) {
                                    gid_sim_start = g_sp_start;
                                    start_pos = {it_ep->second->x_in,
                                                 it_ep->second->y_in, 0.0};
                                    start_yaw = it_ep->second->yaw_in;
                                }
                            }
                        }
                    }
                    target_pt = {ep_last->x_out, ep_last->y_out, 0.0};
                    target_yaw = ep_last->yaw_out;
                }

                // Маска: активны элементы lattice от точки старта симуляции
                // (gid_sim_start = gid_first при изолированном старте, либо
                // первый gid супериода при startFromSpStart=true) до конца окна.
                std::vector<bool> mask(lattice.size(), false);
                for (size_t k = 0; k < meta.size(); ++k) {
                    int gid = meta[k].global_id;
                    if (gid >= gid_sim_start && gid <= gid_mask_hi) mask[k] = true;
                }
                lattice.setActiveMask(mask);

                // Группы текущего окна (из карты). Дедупликация через set --
                // одна и та же multi-SP группа может быть зарегистрирована
                // через несколько gid в одном окне.
                std::set<int> currentWinGroups;
                if (lensWindowing) {
                    // все настраиваемые группы по gid в диапазоне окна (линзы +
                    // настраиваемые магниты между ними, если они в groups).
                    for (int gid : activeGids)
                        if (gid >= gid_first && gid <= gid_last) {
                            auto it = groupsByGid.find(gid);
                            if (it != groupsByGid.end())
                                for (int gi : it->second) currentWinGroups.insert(gi);
                        }
                } else {
                    for (int j = ws; j <= we; ++j) {
                        int gid = activeGids[j];
                        auto it = groupsByGid.find(gid);
                        if (it != groupsByGid.end())
                            for (int gi : it->second) currentWinGroups.insert(gi);
                    }
                }
                // Добавляем в накопительный набор pass'а + фиксируем окно входа.
                // n_new = сколько групп ВПЕРВЫЕ появилось в этом окне.
                int n_new = 0;
                for (int gi : currentWinGroups) {
                    if (!seenGroups.count(gi)) ++n_new;
                    seenGroups.insert(gi);
                    if (!groupEnterWin.count(gi)) groupEnterWin[gi] = win_idx;
                }

                // Активные группы для оптимизации:
                //   accumulateGroups=true  -> seenGroups (все ручки от начала pass'а)
                //   accumulateGroups=false -> только currentWinGroups (стандарт)
                const std::set<int>& chosen = cfg.perElement.accumulateGroups
                                                ? seenGroups : currentWinGroups;
                // lagged-lens: линза (groupLagged) активируется лишь через
                // lensLagWindows окон после первого появления — когда поворотные
                // в её области уже устаканились. Поворотные активны сразу.
                const int lensLag = cfg.perElement.lensLagWindows;
                std::vector<int> activeGroups;
                int lens_frozen = 0;
                activeGroups.reserve(chosen.size());
                for (int gi : chosen) {
                    if (lensLag > 0 && groupLagged[gi]) {
                        auto itE = groupEnterWin.find(gi);
                        int entered = (itE != groupEnterWin.end()) ? itE->second : win_idx;
                        if (win_idx - entered < lensLag) { ++lens_frozen; continue; }
                    }
                    activeGroups.push_back(gi);
                }
                if (activeGroups.empty()) {
                    lattice.clearActiveMask();
                    ws = nextWs(ws, we);
                    continue;  // окно без настраиваемых групп
                }
                // Бюджет окна (единый алгоритм):
                //   total_iters>0   -> РАВНАЯ доля: total_iters / num_windows;
                //   itersPerParam>0 -> itersPerParam × N_active (растёт при accumulate);
                //   иначе           -> фиксированное m_iter.
                int m_iter_dynamic;
                if (budgetPerWindow > 0)
                    m_iter_dynamic = budgetPerWindow;
                else if (cfg.perElement.itersPerParam > 0)
                    m_iter_dynamic = cfg.perElement.itersPerParam * (int)activeGroups.size();
                else
                    m_iter_dynamic = m_iter;

                // Локальные индексы НОВЫХ групп окна (для фокусированной фазы при
                // itersPerNew): первые phase1_iters итераций greedy выбирает группу
                // ТОЛЬКО среди новых -> новый магнит реально настраивается
                // itersPerNew раз. Затем (остаток) -- по всем активным группам.
                std::vector<int> newGroupLocal;
                int phase1_iters = 0;
                if (cfg.perElement.itersPerNew > 0) {
                    for (size_t k = 0; k < activeGroups.size(); ++k) {
                        auto itE = groupEnterWin.find(activeGroups[k]);
                        if (itE != groupEnterWin.end() && itE->second == win_idx)
                            newGroupLocal.push_back((int)k);
                    }
                    phase1_iters = (int)newGroupLocal.size() * cfg.perElement.itersPerNew;
                }

                // Запоминаем стартовые значения активных групп (для отчёта о
                // дельтах в конце окна).
                std::vector<double> winStartVal(activeGroups.size());
                for (size_t k = 0; k < activeGroups.size(); ++k)
                    winStartVal[k] = rgroups[activeGroups[k]].currentValue;

                // Subset magnetPos:
                //   "composite":      магниты с gid в окне + все extra TARGET'ы.
                //   "composite_full": все lattice-магниты с gid ≤ gid_last
                //                     (от старта до конца окна). Extra TARGET'ы
                //                     не включаются -- они после конца окна и
                //                     до них pucok не доходит.
                const bool useComposite = (cfg.perElement.lossMode == "composite");
                const bool useFullTraj  = (cfg.perElement.lossMode == "composite_full");
                std::vector<int> magIdxsInWin;
                std::vector<int> magIdxsFullToWin;
                if (useComposite) {
                    magIdxsInWin.reserve(magnetPos.size());
                    for (size_t mi = 0; mi < magnetPos.size(); ++mi) {
                        int g = (mi < magnetGid.size()) ? magnetGid[mi] : -1;
                        if (g == -1) { magIdxsInWin.push_back((int)mi); continue; }
                        if (g >= gid_first && g <= gid_last)
                            magIdxsInWin.push_back((int)mi);
                    }
                } else if (useFullTraj) {
                    // Все lattice-магниты с gid ≤ gid_last + все extra TARGETs
                    // (gid=-1). Финиш-плоскость теперь глобальная, и target'ы
                    // на drift'е достижимы и важны для оптимизатора.
                    magIdxsFullToWin.reserve(magnetPos.size());
                    for (size_t mi = 0; mi < magnetPos.size(); ++mi) {
                        int g = (mi < magnetGid.size()) ? magnetGid[mi] : -1;
                        if (g == -1) { magIdxsFullToWin.push_back((int)mi); continue; }
                        if (g > 0 && g <= gid_last)
                            magIdxsFullToWin.push_back((int)mi);
                    }
                }
                // LENS-режим: loss-магниты окна. Внутренние окна -- магниты с gid в
                // [gid_first, gid_last] (линзы + между ними), без дальних целей.
                // Последнее окно -- все магниты от gid_first до конца + extra
                // TARGET'ы (gid=-1), т.к. финиш глобальный.
                std::vector<int> magIdxsLens;
                if (lensWindowing) {
                    magIdxsLens.reserve(magnetPos.size());
                    for (size_t mi = 0; mi < magnetPos.size(); ++mi) {
                        int g = (mi < magnetGid.size()) ? magnetGid[mi] : -1;
                        if (g == -1) { if (lensLast) magIdxsLens.push_back((int)mi); continue; }
                        if (g >= gid_first && (lensLast || g <= gid_last))
                            magIdxsLens.push_back((int)mi);
                    }
                }

                // composite_next (вариант B): ЛОКАЛЬНАЯ финиш-плоскость у СЛЕДУЮЩЕГО
                // магнита (gid_last+1). Пучок нанизывается окно-за-окном на опорную
                // орбиту. Внутренние окна: цель = центр след. магнита, локальная
                // плоскость, дальних целей НЕТ. ПОСЛЕДНЕЕ окно: глобальная плоскость
                // (полёт до 50м-мишени) + дальние drift/финиш-цели (gid=-1, вес из
                // JSON). Доминирование дальних целей провоцировало "выстрел" одним
                // квадруполем (QF6B за 2-й линейкой DL2B) в обход орбиты -- лечится
                // НИЗКИМ весом extra-целей (вес 1: нить SP1 28 >> дальние 5).
                const bool useNextMagnet = (cfg.perElement.lossMode == "composite_next");
                // Контрольные точки composite_next: ПЛОТНЫЙ опорный путь (pathPoints),
                // если задан; иначе fallback на центры loss-магнитов (magnetPos).
                const bool usePath = useNextMagnet && !cfg.pathPoints.empty();
                std::vector<int> magIdxsNext;
                FieldVector nextTarget    = target_pt;
                double      nextTargetYaw = target_yaw;
                bool        nextLocalFinish = true;
                if (useNextMagnet) {
                    const int gid_next = gid_last + 1;
                    auto itNep = endpointsByGid.find(gid_next);
                    const bool internalWin = (itNep != endpointsByGid.end());
                    if (internalWin) {
                        // внутреннее окно: цель = центр следующего магнита, ЛОКАЛЬНАЯ плоскость.
                        const auto* nep = itNep->second;
                        nextTarget = {0.5*(nep->x_in + nep->x_out),
                                      0.5*(nep->y_in + nep->y_out), 0.0};
                        nextTargetYaw   = nep->yaw_in;
                        nextLocalFinish = true;
                    } else {
                        // последнее окно: финиш-плоскость ОСТАЁТСЯ глобальной (50м за
                        // drift'ом) -- пучок летит весь путь, нанизываясь на точки пути
                        // (вкл. уплотнённый 50м-дрейф). Плоскость задаёт, ГДЕ остановить
                        // трек; дальних точек-целей с плечом 50м нет.
                        nextLocalFinish = false;
                    }
                    // Подмножество контрольных точек: gid <= gid_last (+ след. магнит).
                    if (usePath) {
                        magIdxsNext.reserve(cfg.pathPoints.size());
                        for (size_t mi = 0; mi < cfg.pathPoints.size(); ++mi) {
                            int g = cfg.pathPointGids[mi];
                            if (g > 0 && (g <= gid_last || (internalWin && g == gid_next)))
                                magIdxsNext.push_back((int)mi);
                        }
                    } else {
                        magIdxsNext.reserve(magnetPos.size());
                        for (size_t mi = 0; mi < magnetPos.size(); ++mi) {
                            int g = (mi < magnetGid.size()) ? magnetGid[mi] : -1;
                            if (g > 0 && (g <= gid_last || (internalWin && g == gid_next)))
                                magIdxsNext.push_back((int)mi);
                        }
                    }
                }

                // ============ PREFIX-OPTIMIZATION ============
                // Если start_from_sp_start=true И accumulate_groups=false:
                // делаем prefix-run ОДИН РАЗ перед окном. Snapshot сохраняется
                // и переиспользуется во всех m итерациях. Это корректно потому,
                // что мутации текущего окна не касаются магнитов ДО окна
                // (mask их выключает, и нет ручек на них в activeGroups).
                //
                // При accumulate_groups=true мутируемые группы МОГУТ влиять
                // на магниты до окна (их параметры тоже меняются). Тогда
                // snapshot, посчитанный ОДИН РАЗ, устаревает -- частица в
                // окне получает неверное стартовое состояние. Чтобы избежать
                // этой ошибки, отключаем prefix-кэш: каждая итерация запускает
                // полную симуляцию от старта. Это медленнее, но корректно.
                std::vector<ParticleState> prefixSnap;
                const std::vector<ParticleState>* snapPtr = nullptr;
                const bool prefix_safe = cfg.perElement.startFromSpStart
                                       && !cfg.perElement.accumulateGroups;
                if (prefix_safe) {
                    if (runPrefixSim(start_pos, start_yaw, prefixSnap)) {
                        snapPtr = &prefixSnap;
                    } else {
                        std::cout << "  WARN: prefix-sim не дошла до entry окна "
                                  << gid_first << " (улетел пучок)\n";
                    }
                }
                // После prefix снова ставим mask окна
                lattice.setActiveMask(mask);

                // Унифицированная обёртка вокруг выбора режима loss.
                // outMean (опц.) -- среднее отклонение траекторий, мм (для гейта).
                // При accumulate_groups + start_from_sp_start prefix пересчитывается
                // на КАЖДОЙ итерации (магниты ДО окна могут меняться через
                // accumulate-набор групп -- старый кэш невалиден).
                auto evalWindow = [&](double* outMean = nullptr) -> double {
                    if (outMean) *outMean = std::numeric_limits<double>::infinity();
                    if (lensWindowing) {
                        // LENS: маска окна уже установлена (активны только элементы
                        // окна; для последнего окна -- до конца lattice). Инжектируем
                        // свежий прямой пучок в start_pos, считаем составной loss
                        // (position + parallelism + envelope_decay) локально в окне.
                        return runFullSimComposite(target_pt, target_yaw, magIdxsLens,
                                                   outMean,
                                                   nullptr, nullptr, nullptr, nullptr,
                                                   lensLocalFinish,
                                                   nullptr, nullptr, nullptr,
                                                   &start_pos, start_yaw,
                                                   cfg.perElement.injectSpreadM);
                    }
                    if (useNextMagnet) {
                        // composite_next: полёт от старта SP до ЛОКАЛЬНОЙ плоскости у
                        // следующего магнита (или глобальной для последнего окна).
                        lattice.clearActiveMask();
                        double L = runFullSimComposite(nextTarget, nextTargetYaw,
                                                        magIdxsNext, outMean,
                                                        nullptr, nullptr, nullptr, nullptr,
                                                        nextLocalFinish,
                                                        usePath ? &cfg.pathPoints       : nullptr,
                                                        usePath ? &cfg.pathPointGids    : nullptr,
                                                        usePath ? &cfg.pathPointWeights : nullptr);
                        lattice.setActiveMask(mask);
                        return L;
                    }
                    if (useFullTraj) {
                        // composite_full: полная симуляция от старта до финиш-плоскости
                        // окна. mask должна быть OFF; восстановим после вызова.
                        lattice.clearActiveMask();
                        double L = runFullSimComposite(target_pt, target_yaw,
                                                        magIdxsFullToWin, outMean);
                        lattice.setActiveMask(mask);
                        return L;
                    }
                    const std::vector<ParticleState>* effSnap = snapPtr;
                    std::vector<ParticleState> dynSnap;
                    if (cfg.perElement.startFromSpStart
                        && cfg.perElement.accumulateGroups) {
                        // mask временно снимаем для prefix-run, потом восстановим
                        lattice.clearActiveMask();
                        if (runPrefixSim(start_pos, start_yaw, dynSnap)) {
                            effSnap = &dynSnap;
                        } else {
                            effSnap = nullptr;  // улетел пучок -- спавн на эталоне
                        }
                        lattice.setActiveMask(mask);
                    }
                    if (useComposite)
                        return runWindowSimComposite(start_pos, start_yaw,
                                                      target_pt, target_yaw,
                                                      magIdxsInWin, effSnap, outMean);
                    return runWindowSim(start_pos, start_yaw, target_pt, target_yaw,
                                         effSnap);
                };

                // Baseline для этого окна
                double winBaselineMean = 0.0;
                double winBaseline = evalWindow(&winBaselineMean);
                double curLoss = winBaseline;
                double winBest = curLoss;
                double winBestMean = winBaselineMean;  // mean при лучшем loss
                int accepted = 0;
                int kicks = 0;
                // Счётчик accept'ов по группам: знаем, какая ручка чаще всего
                // работала в этом окне.
                std::vector<int> acc_per_group(activeGroups.size(), 0);

                // Снимок ЛУЧШИХ значений активных групп -- нужен при киках, т.к.
                // кик безусловно принимается и может оставить окно в худшем
                // состоянии, чем winBest. В конце окна восстанавливаем лучшее.
                std::vector<double> winBestVals(activeGroups.size());
                for (size_t k = 0; k < activeGroups.size(); ++k)
                    winBestVals[k] = rgroups[activeGroups[k]].currentValue;

                const int    kickEvery = cfg.kickEvery;
                const double kickScale = cfg.kickScale;

                // Гейт продвижения окна: пока winBestMean > порога, окно
                // продолжает оптимизироваться (сверх обычного бюджета) до
                // защитного потолка hardCap. 0 => гейт выключен.
                const double advThr = cfg.perElement.advanceMeanThresholdMm;
                const bool   gateOn = (advThr > 0.0);
                const int    hardCap = gateOn
                    ? m_iter_dynamic * std::max(1, cfg.perElement.advanceMaxIterMult)
                    : m_iter_dynamic;

                const bool lineSearch = (cfg.perElement.searchMode == "line_search");
                int it = 0;
                if (lineSearch) {
                    // ============ LINE-SEARCH (последовательное приближение) ============
                    // Покоординатный спуск: проходим по группам окна, каждую
                    // настраиваем детерминированным coarse-to-fine 1D-поиском по
                    // всему допустимому диапазону. K точек на раунд, зум вокруг
                    // ЛУЧШЕЙ известной точки (включая incumbent) с уполовиниванием
                    // интервала. n раундов/элемент распределяется из бюджета окна
                    // равномерно; остаток бюджета -> повторные свипы.
                    const int K  = std::max(3, cfg.perElement.lineSearchPoints);
                    const int nG = (int)activeGroups.size();
                    const int budget = hardCap;  // = m_iter_dynamic (гейт здесь не используем)
                    const int roundsPerGroup = budgetRounds(budget, nG, K,
                                                            cfg.perElement.lineSearchMaxRounds,
                                                            cfg.perElement.lineSearchMinRounds);
                    const int sweepCost = std::max(1, nG * roundsPerGroup * K);
                    const int sweeps = std::max(1, budget / sweepCost);

                    auto commitWinBest = [&](double L, double mean) {
                        if (L < winBest) {
                            winBest = L; winBestMean = mean;
                            for (size_t k = 0; k < activeGroups.size(); ++k)
                                winBestVals[k] = rgroups[activeGroups[k]].currentValue;
                        }
                    };
                    // 1D-поиск по одной группе (остальные группы зафиксированы).
                    auto lineSearchGroup = [&](int gli) {
                        auto& rg = rgroups[activeGroups[gli]];
                        const bool isCoeff = (rg.cfg.target.rfind("COEFF", 0) == 0);
                        const double absLo = isCoeff ? -cfg.perElement.coeffRange
                                                     : -cfg.perElement.currentRangeAmp;
                        const double absHi = isCoeff ?  cfg.perElement.coeffRange
                                                     :  cfg.perElement.currentRangeAmp;
                        // Локальный радиус (против railing): поиск в окрестности текущего.
                        const double localR = isCoeff ? cfg.perElement.lineSearchLocalCoeff
                                                      : cfg.perElement.lineSearchLocalCurrent;
                        const double lo = (localR > 0.0) ? std::max(absLo, rg.currentValue - localR) : absLo;
                        const double hi = (localR > 0.0) ? std::min(absHi, rg.currentValue + localR) : absHi;
                        // incumbent (текущее значение элемента, зажатое в диапазон)
                        double bestVal = std::min(std::max(rg.currentValue, lo), hi);
                        applyGroupValue(rg, bestVal);
                        double bMean = std::numeric_limits<double>::infinity();
                        double bestLoss = evalWindow(&bMean);
                        ++it; ++total_iter; commitWinBest(bestLoss, bMean);
                        const double incumbentLoss = bestLoss;
                        // coarse-to-fine раунды
                        double a = lo, b = hi;
                        for (int r = 0; r < roundsPerGroup; ++r) {
                            const double h = (b - a) / (K - 1);
                            for (int k = 0; k < K; ++k) {
                                const double x = a + k * h;
                                applyGroupValue(rg, x);
                                double mean = std::numeric_limits<double>::infinity();
                                const double L = evalWindow(&mean);
                                ++it; ++total_iter;
                                if (L < bestLoss) { bestLoss = L; bestVal = x; bMean = mean; }
                                commitWinBest(L, mean);
                            }
                            // зум вокруг ЛУЧШЕЙ известной точки (incumbent или sample),
                            // уполовинивая интервал; защищает от увода в плохую область.
                            const double newHalf = h * 0.5;
                            a = std::max(lo, bestVal - newHalf);
                            b = std::min(hi, bestVal + newHalf);
                        }
                        applyGroupValue(rg, bestVal);  // зафиксировать лучшее значение группы
                        curLoss = bestLoss;
                        if (bestLoss < incumbentLoss) { accepted++; acc_per_group[gli]++; }
                    };
                    for (int s = 0; s < sweeps && it < hardCap; ++s)
                        for (int gli = 0; gli < nG && it < hardCap; ++gli)
                            lineSearchGroup(gli);
                } else {
                while (it < hardCap) {
                    ++it;
                    // KICK: каждые kickEvery итераций возмущаем ВСЕ активные
                    // группы окна и безусловно принимаем (выход из локального
                    // минимума). winBest/winBestVals/winBestMean сохраняют лучшее.
                    if (kickEvery > 0 && (it % kickEvery == 0)) {
                        for (int gi : activeGroups) {
                            auto& rg = rgroups[gi];
                            double d = randomStep(rg, kickScale);
                            applyGroupValue(rg, rg.currentValue + d);
                        }
                        double stepMean = std::numeric_limits<double>::infinity();
                        curLoss = evalWindow(&stepMean);
                        if (curLoss < winBest) {
                            winBest = curLoss;
                            winBestMean = stepMean;
                            for (size_t k = 0; k < activeGroups.size(); ++k)
                                winBestVals[k] = rgroups[activeGroups[k]].currentValue;
                        }
                        kicks++;
                        total_iter++;
                    } else {
                        // Обычный шаг: одна группа, greedy accept.
                        // Фаза 1 (первые phase1_iters итер.): только НОВЫЕ группы.
                        // Фаза 2 (остаток): все активные группы.
                        int gi_local;
                        if (cfg.perElement.itersPerNew > 0 && it <= phase1_iters
                            && !newGroupLocal.empty())
                            gi_local = newGroupLocal[std::rand() % (int)newGroupLocal.size()];
                        else
                            gi_local = std::rand() % (int)activeGroups.size();
                        int gi = activeGroups[gi_local];
                        auto& rg = rgroups[gi];
                        double oldVal = rg.currentValue;
                        double delta = randomStep(rg, 1.0);
                        applyGroupValue(rg, oldVal + delta);
                        double stepMean = std::numeric_limits<double>::infinity();
                        double newLoss = evalWindow(&stepMean);
                        if (newLoss < curLoss) {
                            curLoss = newLoss;
                            if (newLoss < winBest) {
                                winBest = newLoss;
                                winBestMean = stepMean;
                                for (size_t k = 0; k < activeGroups.size(); ++k)
                                    winBestVals[k] = rgroups[activeGroups[k]].currentValue;
                            }
                            accepted++;
                            acc_per_group[gi_local]++;
                        } else {
                            applyGroupValue(rg, oldVal);
                        }
                        total_iter++;
                    }

                    // Условие выхода: обычный бюджет отработан + (гейт выкл ИЛИ
                    // mean достиг порога). Без гейта это просто it >= m_iter_dynamic.
                    if (it >= m_iter_dynamic) {
                        if (!gateOn) break;
                        if (winBestMean <= advThr) break;
                        // иначе продолжаем до hardCap
                    }
                }
                }  // конец else (random-режим)
                const int iters_used = it;
                const bool gate_met = (!gateOn) || (winBestMean <= advThr);

                // Восстанавливаем ЛУЧШЕЕ найденное состояние окна (при киках
                // последняя итерация могла оставить пучок в худшем состоянии).
                for (size_t k = 0; k < activeGroups.size(); ++k)
                    applyGroupValue(rgroups[activeGroups[k]], winBestVals[k]);
                curLoss = winBest;
                lattice.clearActiveMask();

                // ОСНОВНОЕ сообщение по окну -- одна строка с компактной сводкой.
                std::cout << "  [pass " << pass << "/" << passes
                          << " | win " << win_idx << "/" << total_windows
                          << " gid=" << gid_first << ".." << gid_last
                          << " grps=" << activeGroups.size()
                          << (cfg.perElement.itersPerNew > 0 ? (" new=" + std::to_string(n_new)) : "")
                          << (lensLag > 0 ? (" lensFrz=" + std::to_string(lens_frozen)) : "")
                          << "]"
                          << "  init=" << std::scientific << std::setprecision(3) << winBaseline
                          << "  last=" << winBest
                          << "  pathmean=" << std::fixed << std::setprecision(2) << winBestMean << "mm"
                          << "  (" << accepted << "/" << iters_used
                          << ", k=" << kicks;
                if (gateOn) std::cout << (gate_met ? ", gate:OK" : ", gate:CAP");
                std::cout << ")"
                          << "  t=" << std::fixed << std::setprecision(1)
                          << totalTimer.sec() << "s\n"
                          << std::defaultfloat;

                // ===== ДИАГНОСТИКА ПО ОКНУ: 3 лосса (pos/dir/env) + сумма, для
                //       ОБОИХ референсов (центры железа vs .lte-орбита). Оптимизация
                //       идёт по primary (cfg.positionReference); вторая версия pos
                //       печатается рядом, чтобы выбрать правильную loss-функцию.
                if (useFullTraj || useNextMagnet) {
                    lattice.clearActiveMask();
                    double dMean=0, dPos=0, dDir=0, dEnv=0, dPosAlt=0;
                    if (useNextMagnet)
                        runFullSimComposite(nextTarget, nextTargetYaw, magIdxsNext,
                                            &dMean, &dPos, &dDir, &dEnv, &dPosAlt,
                                            nextLocalFinish,
                                            usePath ? &cfg.pathPoints       : nullptr,
                                            usePath ? &cfg.pathPointGids    : nullptr,
                                            usePath ? &cfg.pathPointWeights : nullptr);
                    else
                        runFullSimComposite(target_pt, target_yaw, magIdxsFullToWin,
                                            &dMean, &dPos, &dDir, &dEnv, &dPosAlt);
                    const double wd = cfg.dirWeight, we = cfg.envWeight;
                    const bool primIsOrbit = (cfg.positionReference == "orbit");
                    const char* primName = primIsOrbit ? ".lte " : "iron ";
                    const char* altName  = primIsOrbit ? "iron " : ".lte ";
                    std::cout << std::scientific << std::setprecision(3)
                      << "         loss  dir=" << dDir << "  env=" << dEnv
                      << "  (wd=" << wd << " we=" << we << ")\n"
                      << "           " << primName << "[opt]: pos=" << dPos
                      << "  sum=" << (dPos    + wd*dDir + we*dEnv) << "\n"
                      << "           " << altName  << "      : pos=" << dPosAlt
                      << "  sum=" << (dPosAlt + wd*dDir + we*dEnv) << "\n"
                      << std::defaultfloat;
                }
                if (lensWindowing) {
                    // Диагностика по окну: позиция vs (параллельность+огибающая).
                    lattice.setActiveMask(mask);
                    double dMean=0, dPos=0, dDir=0, dEnv=0;
                    double tot = runFullSimComposite(target_pt, target_yaw, magIdxsLens,
                                        &dMean, &dPos, &dDir, &dEnv, nullptr,
                                        lensLocalFinish, nullptr, nullptr, nullptr,
                                        &start_pos, start_yaw, cfg.perElement.injectSpreadM);
                    lattice.clearActiveMask();
                    std::cout << std::scientific << std::setprecision(3)
                      << "         pos=" << dPos
                      << "  parallel+decay=" << (tot - dPos)
                      << "  (всего=" << tot << ")\n" << std::defaultfloat;
                }
                // Запись в общий history (для combined-режима и пост-анализа).
                // Одна запись на окно с финальным win_best loss.
                history.push_back(HistEntry{
                    (int)total_iter, "PE_WIN", 0.0, 0.0,
                    winBest, winBest, 0.0, 0.0,
                    std::string("WIN_gid") + std::to_string(gid_first) +
                    "_" + std::to_string(gid_last)});

                ws = nextWs(ws, we);   // продвижение окна (пропуск нутра дип. линейки)
            }

            // Запись в history: метка завершения пасса.
            history.push_back(HistEntry{
                (int)total_iter, "PE_PASS", 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0,
                std::string("PASS_END_") + std::to_string(pass)});
        }

        // ============ ФАЗА РАСШИРЕНИЯ В ДРЕЙФ (прогрессивная) ============
        // После укладки нити SP1 окнами -- "вползаем" в дрейф по driftChunkM метров.
        // Рубеж растёт; на каждом контрольные точки = весь путь ДО рубежа (SP1 +
        // дрейф), финиш-плоскость на рубеже, настраиваются ВСЕ группы (line-search,
        // локальный радиус). Пучок обязан следовать линии непрерывно до рубежа ->
        // ни "выстрела", ни расщепления (между точками теперь контроль).
        if (cfg.perElement.driftExtensionM > 0.0 && !cfg.pathPoints.empty()) {
            const int maxLatGid = activeGids.empty() ? 0 : activeGids.back();
            // Виртуальные gid'ы дрейфа -> рубеж-позиция (самая дальняя точка с этим gid).
            std::map<int, FieldVector> driftExtent;
            for (size_t mi = 0; mi < cfg.pathPoints.size(); ++mi) {
                int g = cfg.pathPointGids[mi];
                if (g > maxLatGid) {
                    auto it = driftExtent.find(g);
                    if (it == driftExtent.end() || cfg.pathPoints[mi].x > it->second.x)
                        driftExtent[g] = cfg.pathPoints[mi];
                }
            }
            const double driftYaw = std::atan2(cfg.dirReference.y, cfg.dirReference.x);
            const int    nG = (int)rgroups.size();
            const int    K  = std::max(3, cfg.perElement.lineSearchPoints);
            const int rounds = budgetRounds(cfg.perElement.driftItersPerChunk, nG, K,
                                            cfg.perElement.lineSearchMaxRounds,
                                            cfg.perElement.lineSearchMinRounds);
            std::cout << "\n[DRIFT EXTENSION] чанков=" << driftExtent.size()
                      << "  шаг=" << cfg.perElement.driftChunkM << "м до "
                      << cfg.perElement.driftExtensionM << "м  rounds/grp=" << rounds << "\n";
            // Оценка loss на рубеже: полёт от старта SP до локальной плоскости рубежа,
            // контрольные точки = path с gid<=gidLimit.
            auto evalDrift = [&](const FieldVector& extentPos, const std::vector<int>& mag,
                                 double* outMean) -> double {
                lattice.clearActiveMask();
                return runFullSimComposite(extentPos, driftYaw, mag, outMean,
                                           nullptr, nullptr, nullptr, nullptr, true,
                                           &cfg.pathPoints, &cfg.pathPointGids,
                                           &cfg.pathPointWeights);
            };
            int chunkNo = 0;
            for (auto& kv : driftExtent) {
                const int gidLimit = kv.first;
                const FieldVector extentPos = kv.second;
                ++chunkNo;
                std::vector<int> mag;
                mag.reserve(cfg.pathPoints.size());
                for (size_t mi = 0; mi < cfg.pathPoints.size(); ++mi)
                    if (cfg.pathPointGids[mi] > 0 && cfg.pathPointGids[mi] <= gidLimit)
                        mag.push_back((int)mi);
                // line-search по всем группам (1 свип): локальный радиус вокруг текущего.
                for (int gi = 0; gi < nG; ++gi) {
                    auto& rg = rgroups[gi];
                    const bool isCoeff = (rg.cfg.target.rfind("COEFF", 0) == 0);
                    const double absLo = isCoeff ? -cfg.perElement.coeffRange : -cfg.perElement.currentRangeAmp;
                    const double absHi = isCoeff ?  cfg.perElement.coeffRange :  cfg.perElement.currentRangeAmp;
                    const double localR = isCoeff ? cfg.perElement.lineSearchLocalCoeff
                                                  : cfg.perElement.lineSearchLocalCurrent;
                    const double lo = (localR > 0.0) ? std::max(absLo, rg.currentValue - localR) : absLo;
                    const double hi = (localR > 0.0) ? std::min(absHi, rg.currentValue + localR) : absHi;
                    double bestVal = std::min(std::max(rg.currentValue, lo), hi);
                    applyGroupValue(rg, bestVal);
                    double m0 = 0.0;
                    double bestLoss = evalDrift(extentPos, mag, &m0);
                    ++total_iter;
                    double a = lo, b = hi;
                    for (int r = 0; r < rounds; ++r) {
                        const double h = (b - a) / (K - 1);
                        for (int k = 0; k < K; ++k) {
                            const double x = a + k * h;
                            applyGroupValue(rg, x);
                            double mm = 0.0;
                            const double L = evalDrift(extentPos, mag, &mm);
                            ++total_iter;
                            if (L < bestLoss) { bestLoss = L; bestVal = x; }
                        }
                        const double nh = h * 0.5;
                        a = std::max(lo, bestVal - nh);
                        b = std::min(hi, bestVal + nh);
                    }
                    applyGroupValue(rg, bestVal);
                }
                double fMean = 0.0;
                double fL = evalDrift(extentPos, mag, &fMean);
                if (chunkNo == 1 || chunkNo % 5 == 0 || gidLimit == driftExtent.rbegin()->first)
                    std::cout << "  [drift " << chunkNo << "/" << driftExtent.size()
                              << " gid<=" << gidLimit << "]  loss=" << std::scientific
                              << std::setprecision(3) << fL << "  pathmean=" << std::fixed
                              << std::setprecision(2) << fMean << "mm  t=" << std::setprecision(1)
                              << totalTimer.sec() << "s\n" << std::defaultfloat;
            }
            lattice.clearActiveMask();
        }

        // ---- Финальная глобальная оценка ПОСЛЕ всех passes ----
        {
            lattice.clearActiveMask();
            double pgL, pgMean, pgMax, pos_c, dir_c, env_c;
            bool after_ok = runTrackingAndLoss(pgL, pgMean, pgMax, &pos_c, &dir_c, &env_c);
            if (after_ok) {
                bestLoss = pgL;   // реальный финальный loss оконного режима (для хвостового "Best loss")
                std::cout << "\n[global FINAL]  composite=" << std::scientific
                          << std::setprecision(3) << pgL
                          << "  pos=" << pos_c << "  dir=" << dir_c << "  env=" << env_c
                          << "  mean=" << std::fixed << std::setprecision(2)
                          << pgMean << "mm  max=" << pgMax << "mm  t="
                          << std::setprecision(1) << totalTimer.sec() << "s\n"
                          << std::defaultfloat;
            } else {
                std::cout << "\n[global FINAL]  -- симуляция не дошла до финиша\n";
            }

            // ============== РЕЗЮМЕ: before vs after ==============
            if (before_ok && after_ok) {
                auto pct = [](double a, double b) -> double {
                    return (std::abs(a) > 1e-30) ? 100.0 * (a - b) / a : 0.0;
                };
                std::cout << "\n========== ИТОГ (оконная оптимизация) ==========\n";
                // Метод / loss / путь / дрейф / бюджет.
                std::cout << "  метод: "
                          << (cfg.perElement.searchMode == "line_search"
                                  ? "line-search (последовательное приближение)"
                                  : "random (свипы)");
                if (cfg.perElement.lineSearchLocalCurrent > 0.0)
                    std::cout << ", локальный радиус ±" << std::setprecision(0)
                              << std::fixed << cfg.perElement.lineSearchLocalCurrent << "A";
                if (cfg.kickEvery > 0) std::cout << ", кики каждые " << cfg.kickEvery;
                std::cout << "\n  loss: " << cfg.perElement.lossMode;
                if (!cfg.pathPoints.empty())
                    std::cout << " (опорный путь " << cfg.pathPoints.size() << " точек)";
                if (cfg.perElement.driftExtensionM > 0.0)
                    std::cout << " + вползание в дрейф до " << std::setprecision(0)
                              << std::fixed << cfg.perElement.driftExtensionM << "м";
                std::cout << "\n  бюджет: " << total_iter << " симуляций, "
                          << passes << " проход(ов)\n";
                std::cout << std::scientific << std::setprecision(3);
                std::cout << "  composite:  " << before_L   << "  ->  " << pgL
                          << "    (" << std::fixed << std::setprecision(2)
                          << pct(before_L, pgL) << "%)\n";
                std::cout << std::scientific << std::setprecision(3);
                std::cout << "  position:   " << before_pos << "  ->  " << pos_c
                          << "    (" << std::fixed << std::setprecision(2)
                          << pct(before_pos, pos_c) << "%)\n";
                std::cout << std::scientific << std::setprecision(3);
                std::cout << "  direction:  " << before_dir << "  ->  " << dir_c
                          << "    (" << std::fixed << std::setprecision(2)
                          << pct(before_dir, dir_c) << "%)\n";
                std::cout << std::scientific << std::setprecision(3);
                std::cout << "  envelope:   " << before_env << "  ->  " << env_c
                          << "    (" << std::fixed << std::setprecision(2)
                          << pct(before_env, env_c) << "%)\n";
                std::cout << std::fixed << std::setprecision(3);
                std::cout << "  mean:       " << before_mean << " mm  ->  "
                          << pgMean << " mm    (" << std::setprecision(2)
                          << pct(before_mean, pgMean) << "%)\n";
                std::cout << std::fixed << std::setprecision(3);
                std::cout << "  max:        " << before_max  << " mm  ->  "
                          << pgMax << " mm    (" << std::setprecision(2)
                          << pct(before_max, pgMax) << "%)\n";
                // Сводка по токам (физичность решения): медиана/макс/на упоре.
                std::vector<double> absI;
                for (const auto& rg : rgroups)
                    if (rg.cfg.target.rfind("COEFF", 0) != 0 &&
                        std::abs(rg.currentValue) > 1e-9)
                        absI.push_back(std::abs(rg.currentValue));
                if (!absI.empty()) {
                    std::sort(absI.begin(), absI.end());
                    int railed = 0;
                    for (double v : absI)
                        if (v >= cfg.perElement.currentRangeAmp - 0.1) ++railed;
                    std::cout << std::fixed << std::setprecision(1)
                              << "  токи:       " << absI.size() << " ненул., медиана|I|="
                              << absI[absI.size()/2] << "A, max|I|=" << absI.back()
                              << "A, на упоре=" << railed << "\n";
                }
                std::cout << "  время:      " << std::fixed << std::setprecision(1)
                          << totalTimer.sec() << "s\n" << std::defaultfloat;
                std::cout << "================================================\n";
            }
        }
        lattice.clearActiveMask();

        // Сохраняем лучшие значения и пропускаем MC/gradient ниже,
        // переходя к финализации (запись sp_elements.csv, history, trajectory).
        for (size_t k = 0; k < rgroups.size(); ++k) bestValues[k] = rgroups[k].currentValue;
        goto skip_mc_loop;
    }

    // ========================================================================
    // GRADIENT MODE — численный градиент через конечные разности + Adam +
    // backtracking line search. Делит интерфейс с MC: те же applyGroupValue,
    // runTrackingAndLoss, history, bestLoss/bestValues, writeHistory.
    // ========================================================================
    if (cfg.method == "gradient") {
        const auto& gc = cfg.gradient;
        std::cout << "\n[GRADIENT MODE] params=" << rgroups.size()
                  << "  fd=" << gc.fdMethod
                  << "  rule=" << gc.updateRule
                  << "  lr=" << gc.learningRate;
        if (gc.updateRule == "adam") {
            std::cout << "  β1=" << gc.beta1 << " β2=" << gc.beta2;
        }
        std::cout << "  ε_abs=" << gc.fdStepAbs << " ε_pct=" << gc.fdStepPercent
                  << "  max_halvings=" << gc.maxHalvings << "\n";

        // ε_j по типу группы: для arg_val — fdStepAbs (А), для COEFF — fdStepPercent · |θ_j|.
        auto fdEps = [&](const RuntimeGroup& rg) -> double {
            if (rg.cfg.target == "COEFF_HW" || rg.cfg.target == "COEFF_INST") {
                double eps = gc.fdStepPercent * std::abs(rg.currentValue);
                if (eps < 1e-12) eps = gc.fdStepPercent;  // защита от θ=0
                return eps;
            }
            return gc.fdStepAbs;
        };

        const size_t M = rgroups.size();
        std::vector<double> mAdam(M, 0.0), vAdam(M, 0.0);
        std::vector<double> grad(M, 0.0);
        std::vector<double> theta(M, 0.0);
        for (size_t j = 0; j < M; ++j) theta[j] = rgroups[j].currentValue;

        bool central = (gc.fdMethod == "central");
        int totalFdCalls = 0, totalLsCalls = 0;

        for (int it = 1; it <= gc.nIter; ++it) {
            StepTimer itTimer;
            int failedAxes = 0;
            int fdCalls = 0;

            // ---- A. Численный градиент ----
            // Базовое значение L_base нужно только для forward FD; для central не используется.
            double L_base = curLoss;
            for (size_t j = 0; j < M; ++j) {
                double eps_j = fdEps(rgroups[j]);
                double theta_j = theta[j];

                bool fail = false;
                double L_plus = 0, L_minus = 0;
                double dummyM = 0, dummyX = 0;

                applyGroupValue(rgroups[j], theta_j + eps_j);
                if (!runTrackingAndLoss(L_plus, dummyM, dummyX)) fail = true;
                fdCalls++;

                if (central) {
                    applyGroupValue(rgroups[j], theta_j - eps_j);
                    if (!runTrackingAndLoss(L_minus, dummyM, dummyX)) fail = true;
                    fdCalls++;
                }
                // Восстанавливаем θ_j (важно: до этого мы сместили его внутри FD).
                applyGroupValue(rgroups[j], theta_j);

                if (fail) {
                    grad[j] = 0.0;             // не двигаемся по оси с failure
                    failedAxes++;
                } else if (central) {
                    grad[j] = (L_plus - L_minus) / (2.0 * eps_j);
                } else {
                    grad[j] = (L_plus - L_base) / eps_j;
                }
            }
            totalFdCalls += fdCalls;

            // Норма градиента + индекс максимальной по модулю компоненты (для диагностики)
            double gNorm2 = 0, maxAbsG = 0;
            int    maxIdx = -1;
            for (size_t j = 0; j < M; ++j) {
                gNorm2 += grad[j]*grad[j];
                double a = std::abs(grad[j]);
                if (a > maxAbsG) { maxAbsG = a; maxIdx = (int)j; }
            }
            double gNorm = std::sqrt(gNorm2);

            // ---- B. Update step: normalized_gd | adam ----
            std::vector<double> dTheta(M, 0.0);
            if (gc.updateRule == "normalized_gd") {
                // Δθ_j = -lr · g_j / max_k|g_k|.
                // Параметр с максимальным |g| смещается на ±lr, остальные пропорционально.
                if (maxAbsG > 1e-20) {
                    for (size_t j = 0; j < M; ++j)
                        dTheta[j] = -gc.learningRate * grad[j] / maxAbsG;
                }
                // Adam-моменты тоже обновляем — на случай если потом переключатся на adam:
                // нет, не делаем. Adam-state живёт изолированно.
            } else if (gc.updateRule == "adam") {
                // Bias-correction: m̂ = m / (1-β1^t), v̂ = v / (1-β2^t).
                double biasCorr1 = 1.0 - std::pow(gc.beta1, (double)it);
                double biasCorr2 = 1.0 - std::pow(gc.beta2, (double)it);
                for (size_t j = 0; j < M; ++j) {
                    mAdam[j] = gc.beta1 * mAdam[j] + (1 - gc.beta1) * grad[j];
                    vAdam[j] = gc.beta2 * vAdam[j] + (1 - gc.beta2) * grad[j] * grad[j];
                    double mHat = mAdam[j] / biasCorr1;
                    double vHat = vAdam[j] / biasCorr2;
                    dTheta[j] = -gc.learningRate * mHat / (std::sqrt(vHat) + gc.epsAdam);
                }
            } else {
                std::cerr << "ERROR: unknown gradient.update_rule '"
                          << gc.updateRule << "'  (expected 'normalized_gd' or 'adam')\n";
                return 1;
            }

            // ---- C. Backtracking line search ----
            double L_old = curLoss;
            double L_new = L_old;
            double meanMM_new = curMeanMM, maxMM_new = curMaxMM;
            double alpha = 1.0;
            double alphaUsed = 0.0;
            double L_try_last = L_old;     // последнее наблюдённое значение (или L_old если все FAIL)
            bool accepted_step = false;
            int lsCalls = 0, lsFails = 0;
            for (int k = 0; k <= gc.maxHalvings; ++k) {
                // θ_new = θ + α · Δθ
                for (size_t j = 0; j < M; ++j) {
                    applyGroupValue(rgroups[j], theta[j] + alpha * dTheta[j]);
                }
                double L_try = 0, mm = 0, mx = 0;
                bool ok = runTrackingAndLoss(L_try, mm, mx);
                lsCalls++;
                alphaUsed = alpha;
                if (!ok) {
                    lsFails++;
                } else {
                    L_try_last = L_try;
                    if (L_try < L_old) {
                        L_new = L_try;
                        meanMM_new = mm;  maxMM_new = mx;
                        for (size_t j = 0; j < M; ++j) theta[j] += alpha * dTheta[j];
                        accepted_step = true;
                        break;
                    }
                }
                alpha *= 0.5;
            }
            if (!accepted_step) L_new = L_try_last;   // покажем последнюю достигнутую L, не stale L_old
            totalLsCalls += lsCalls;

            if (!accepted_step) {
                // Все попытки в line search провалились — откатываем θ к началу итерации.
                for (size_t j = 0; j < M; ++j) applyGroupValue(rgroups[j], theta[j]);
                rejected++;
            } else {
                curLoss = L_new;  curMeanMM = meanMM_new;  curMaxMM = maxMM_new;
                accepted++;
                if (L_new < bestLoss) {
                    bestLoss = L_new;
                    for (size_t j = 0; j < M; ++j) bestValues[j] = theta[j];
                    kicksWithoutBest = 0;
                }
            }

            // ---- D. ЛОГ ВСЕХ ИТЕРАЦИЙ (как просил пользователь) ----
            const char* status = accepted_step
                ? (L_new < bestLoss + 1e-30 && bestLoss == L_new ? "ACCEPT+BEST" : "ACCEPT")
                : "REJECT";
            double dL = L_new - L_old;
            // Имя и знак доминирующей по |g| компоненты — диагностика направления.
            std::string topName = (maxIdx >= 0 && maxIdx < (int)rgroups.size())
                                ? rgroups[maxIdx].cfg.name : "?";
            char topSign = (maxIdx >= 0 && grad[maxIdx] > 0) ? '+' : '-';

            std::cout << "[g " << it << "/" << gc.nIter << "] "
                      << "|g|=" << std::scientific << std::setprecision(3) << gNorm
                      << "  L=" << std::setprecision(6) << L_old
                      << " → " << L_new
                      << "  ΔL=" << std::showpos << std::setprecision(3) << dL << std::noshowpos
                      << "  α=" << std::fixed << std::setprecision(4) << alphaUsed
                      << "  best=" << std::scientific << std::setprecision(6) << bestLoss
                      << "  fd=" << fdCalls << " ls=" << lsCalls
                      << "  fd_fail=" << failedAxes
                      << " ls_fail=" << lsFails
                      << "  top=" << topSign << topName
                      << "  t=" << std::fixed << std::setprecision(1) << itTimer.sec() << "s"
                      << "  [" << status << "]"
                      << std::resetiosflags(std::ios_base::floatfield) << "\n";

            // История + чекпойнт. mark кодируем: GRAD_ACCEPT / GRAD_REJECT / GRAD_BEST.
            std::string mark = accepted_step
                ? (L_new == bestLoss ? "GRAD_BEST" : "GRAD_ACCEPT")
                : "GRAD_REJECT";
            HistEntry e{it, std::string("GRAD"), L_old, L_new, L_new, bestLoss,
                        meanMM_new, maxMM_new, mark};
            history.push_back(e);

            if (it % CHECKPOINT_EVERY == 0) writeHistory(false);
        }

        std::cout << "\n[GRADIENT MODE] done."
                  << "  total FD calls=" << totalFdCalls
                  << "  total LS calls=" << totalLsCalls
                  << "  accepted=" << accepted
                  << "  rejected=" << rejected << "\n";

        // Пропускаем MC-цикл ниже: переход прямо к восстановлению best state.
        goto skip_mc_loop;
    }

    // ============== MAIN MC LOOP (wrapped in {} to prevent goto-crossing init) ==============
    {
    // ============== DYNAMIC LOSS init ==============
    // Период T для синусоидальной модуляции весов composite-loss. Если не задан --
    // авто = 5 × N_params (каждый параметр в среднем мутирует ~5 раз за оборот).
    int dynT = cfg.dynamicLoss.periodIterations;
    if (cfg.dynamicLoss.enabled && dynT <= 0) {
        dynT = std::max(1, 5 * (int)rgroups.size());
    }
    if (cfg.dynamicLoss.enabled) {
        std::cout << "[DYNAMIC LOSS] enabled  amplitude=" << cfg.dynamicLoss.amplitude
                  << "  period=" << dynT << " iters (phase rotates pos/dir/env)\n";
    }
    // Лямбда: для итерации it возвращает три весовых фактора w_pos, w_dir, w_env.
    auto dynamicWeights = [&](int it, double& wp, double& wd, double& we) {
        if (!cfg.dynamicLoss.enabled || dynT <= 0) {
            wp = wd = we = 1.0;
            return;
        }
        const double TWO_PI = 6.283185307179586;
        double th = TWO_PI * (double)it / (double)dynT;
        double A  = cfg.dynamicLoss.amplitude;
        wp = 1.0 + A * std::cos(th);                      // phase 0
        wd = 1.0 + A * std::cos(th - TWO_PI / 3.0);       // phase 2π/3
        we = 1.0 + A * std::cos(th - 2.0 * TWO_PI / 3.0); // phase 4π/3
        // Защита от негативных весов (если A > 1)
        if (wp < 0.0) wp = 0.0;
        if (wd < 0.0) wd = 0.0;
        if (we < 0.0) we = 0.0;
    };
    // Лямбда: комбинированный loss = w_pos·pos + w_dir·dirWeight·dir + w_env·envWeight·env.
    auto composeLoss = [&](double wp, double wd, double we,
                            double pos_c, double dir_c, double env_c) -> double {
        return wp * pos_c
             + wd * cfg.dirWeight * dir_c
             + we * cfg.envWeight * env_c;
    };

    const bool svdGlobal = (cfg.globalSearchMode == "svd");
    const bool lsGlobal  = (cfg.globalSearchMode == "line_search");
    if (svdGlobal) {
        // ====== RESPONSE MATRIX + SVD (Gauss-Newton orbit correction) ======
        // Невязка r = sqrt(w_m)·d(i,m) по всем (частица,магнит). Якобиан J = ∂r/∂param
        // строится конечными разностями (nG симуляций). Шаг Гаусса-Ньютона решает
        // J·Δp ≈ −r усечённым SVD (через собств. разложение нормальной матрицы JᵀJ:
        // отбрасываем малые сингулярные числа = плохо обусловленные комбинации токов).
        // Демпфированный шаг с backtracking по РЕАЛЬНОМУ composite-loss; границы ОДЗ.
        const int nG = (int)rgroups.size();
        const int budget = std::max(1, cfg.perElement.totalIters > 0
                                          ? cfg.perElement.totalIters : cfg.nIter);
        auto isCoeffG = [](const RuntimeGroup& rg){ return rg.cfg.target.rfind("COEFF",0)==0; };
        std::vector<double> du(nG), lo(nG), hi(nG);
        for (int j = 0; j < nG; ++j) {
            // Шаг конечной разности крупный (выше ~0.2мм дискретизац. джиттера closest-approach):
            // ток 5A, coeff 0.05. Слишком мелкий шаг -> якобиан тонет в шуме метрики.
            if (isCoeffG(rgroups[j])) { du[j]=0.05; lo[j]=-cfg.perElement.coeffRange;       hi[j]=cfg.perElement.coeffRange; }
            else                      { du[j]=5.0;  lo[j]=-cfg.perElement.currentRangeAmp;  hi[j]=cfg.perElement.currentRangeAmp; }
        }
        auto clampG = [&](int j, double v){ return std::min(hi[j], std::max(lo[j], v)); };

        // Симметричное собств. разложение (циклический Якоби). A (n×n row-major,
        // перезаписывается); w[n] — собств.числа (=σ²); V (n×n, столбцы — собств.векторы).
        auto jacobiEigen = [](std::vector<double>& A, int n,
                              std::vector<double>& w, std::vector<double>& V) {
            V.assign((size_t)n*n, 0.0);
            for (int i = 0; i < n; ++i) V[(size_t)i*n+i] = 1.0;
            for (int sweep = 0; sweep < 100; ++sweep) {
                double off = 0.0;
                for (int p = 0; p < n; ++p)
                    for (int q = p+1; q < n; ++q) off += A[(size_t)p*n+q]*A[(size_t)p*n+q];
                if (off < 1e-30) break;
                for (int p = 0; p < n; ++p) for (int q = p+1; q < n; ++q) {
                    double apq = A[(size_t)p*n+q];
                    if (std::abs(apq) < 1e-300) continue;
                    double app = A[(size_t)p*n+p], aqq = A[(size_t)q*n+q];
                    double phi = 0.5 * std::atan2(2.0*apq, app - aqq);
                    double c = std::cos(phi), s = std::sin(phi);
                    for (int k = 0; k < n; ++k) {                 // A <- A·G (столбцы p,q)
                        double akp = A[(size_t)k*n+p], akq = A[(size_t)k*n+q];
                        A[(size_t)k*n+p] = c*akp - s*akq;
                        A[(size_t)k*n+q] = s*akp + c*akq;
                    }
                    for (int k = 0; k < n; ++k) {                 // A <- Gᵀ·A (строки p,q)
                        double apk = A[(size_t)p*n+k], aqk = A[(size_t)q*n+k];
                        A[(size_t)p*n+k] = c*apk - s*aqk;
                        A[(size_t)q*n+k] = s*apk + c*aqk;
                    }
                    for (int k = 0; k < n; ++k) {                 // V <- V·G
                        double vkp = V[(size_t)k*n+p], vkq = V[(size_t)k*n+q];
                        V[(size_t)k*n+p] = c*vkp - s*vkq;
                        V[(size_t)k*n+q] = s*vkp + c*vkq;
                    }
                }
            }
            w.assign(n, 0.0);
            for (int i = 0; i < n; ++i) w[i] = A[(size_t)i*n+i];
        };

        std::cout << "[SVD / RESPONSE-MATRIX] groups=" << nG
                  << "  бюджет=" << budget << " симуляций"
                  << "  ОДЗ: |ток|<=" << cfg.perElement.currentRangeAmp
                  << "A, |coeff|<=" << cfg.perElement.coeffRange
                  << (cfg.svdAngleWeight > 0.0
                        ? "  +ПЕРИОДИЧНОСТЬ (выход∥орбите, вес "
                            + std::to_string(cfg.svdAngleWeight) + ")"
                        : "  (только позиция)")
                  << "\n";

        std::vector<double> r0, rj;
        int evalNo = 0;
        const double gBeforeMean = curMeanMM, gBeforeMax = curMaxMM;
        for (int iter = 1; ; ++iter) {
            double L0, m0, x0;
            if (!runTrackingAndLoss(L0, m0, x0, nullptr,nullptr,nullptr,nullptr,nullptr, &r0)) {
                std::cout << "  [SVD] базовый трек не дошёл до финиша — стоп\n"; break;
            }
            ++evalNo;
            const int M = (int)r0.size();
            double rn0 = 0.0; for (double v : r0) rn0 += v*v;   // ||r||² — её и минимизируем
            // Якобиан J (M×nG), по столбцам, конечными разностями.
            std::vector<std::vector<double>> Jc(nG);
            for (int j = 0; j < nG; ++j) {
                double basev = rgroups[j].currentValue;
                double pv = clampG(j, basev + du[j]);
                if (pv == basev) pv = clampG(j, basev - du[j]);
                double step = pv - basev;
                applyGroupValue(rgroups[j], pv);
                double Lj, mj, xj;
                bool ok = runTrackingAndLoss(Lj, mj, xj, nullptr,nullptr,nullptr,nullptr,nullptr, &rj);
                applyGroupValue(rgroups[j], basev);
                ++evalNo;
                Jc[j].assign(M, 0.0);
                if (ok && std::abs(step) > 1e-30 && (int)rj.size() == M)
                    for (int k = 0; k < M; ++k) Jc[j][k] = (rj[k] - r0[k]) / step;
            }
            // A = JᵀJ (nG×nG), g = Jᵀr0 (nG).
            std::vector<double> A((size_t)nG*nG, 0.0), g(nG, 0.0);
            for (int a = 0; a < nG; ++a) {
                const auto& Ja = Jc[a];
                double ga = 0.0; for (int k = 0; k < M; ++k) ga += Ja[k]*r0[k];
                g[a] = ga;
                for (int b = a; b < nG; ++b) {
                    const auto& Jb = Jc[b];
                    double s = 0.0; for (int k = 0; k < M; ++k) s += Ja[k]*Jb[k];
                    A[(size_t)a*nG+b] = s; A[(size_t)b*nG+a] = s;
                }
            }
            std::vector<double> w, V;
            jacobiEigen(A, nG, w, V);
            double wmax = 0.0; for (double v : w) wmax = std::max(wmax, v);
            if (wmax <= 0.0) { std::cout << "  [SVD] вырожденный якобиан — стоп\n"; break; }
            // SVD-УСЕЧЕНИЕ: оставляем только хорошо обусловленные моды (σ²/σ²_max > relTol).
            // Малые сингулярные числа = почти-вырожденные комбинации токов, дающие огромный
            // мусорный шаг 1/w_k — их отбрасываем (суть коррекции орбиты по SVD).
            const double relTol = 1e-3;
            const double tol = relTol * wmax;
            std::vector<double> dp(nG, 0.0);
            int usedModes = 0;
            for (int k = 0; k < nG; ++k) {
                if (w[k] <= tol) continue;
                double vg = 0.0;
                for (int i = 0; i < nG; ++i) vg += V[(size_t)i*nG+k]*g[i];
                double coef = -vg / w[k];
                for (int i = 0; i < nG; ++i) dp[i] += coef * V[(size_t)i*nG+k];
                ++usedModes;
            }
            // ДОВЕРИТЕЛЬНАЯ ОБЛАСТЬ: масштабируем весь шаг так, чтобы за одну итерацию
            // |Δток|<=20A, |Δcoeff|<=0.5 (против переоценки линеаризации на большом шаге).
            double scale = 1.0;
            for (int j = 0; j < nG; ++j) {
                double mx = isCoeffG(rgroups[j]) ? 0.5 : 20.0;
                double ratio = std::abs(dp[j]) / mx;
                if (ratio > scale) scale = ratio;
            }
            if (scale > 1.0) for (int j = 0; j < nG; ++j) dp[j] /= scale;
            // Backtracking вдоль Δp по НЕВЯЗКЕ ||r||² (позиция + периодичность) — её и
            // минимизирует Гаусс-Ньютон. Демпфирование уполовиниванием α.
            std::vector<double> basep(nG);
            for (int j = 0; j < nG; ++j) basep[j] = rgroups[j].currentValue;
            double alpha = 1.0; bool improved = false;
            double stepLoss = L0, stepMean = m0, stepMax = x0, stepRn = rn0;
            for (int bt = 0; bt < 8 && evalNo < budget + nG; ++bt, alpha *= 0.5) {
                for (int j = 0; j < nG; ++j) applyGroupValue(rgroups[j], clampG(j, basep[j] + alpha*dp[j]));
                double Lt, mt, xt;
                bool ok = runTrackingAndLoss(Lt, mt, xt, nullptr,nullptr,nullptr,nullptr,nullptr, &rj);
                ++evalNo;
                double rnt = 0.0; if (ok) for (double v : rj) rnt += v*v;
                if (ok && rnt < rn0) {
                    stepLoss = Lt; stepMean = mt; stepMax = xt; stepRn = rnt; improved = true; break;
                }
            }
            if (improved) {
                curLoss = stepLoss; curMeanMM = stepMean; curMaxMM = stepMax;
                // ||r||² упал => принимаем как лучшее SVD-состояние (эпилог применит bestValues).
                bestLoss = curLoss;
                for (int j = 0; j < nG; ++j) bestValues[j] = rgroups[j].currentValue;
            } else {
                for (int j = 0; j < nG; ++j) applyGroupValue(rgroups[j], basep[j]);  // откат
            }
            std::cout << "[SVD итер " << iter << "]  модов=" << usedModes << "/" << nG
                      << "  α=" << std::fixed << std::setprecision(3) << alpha
                      << "  |r|=" << std::scientific << std::setprecision(3) << std::sqrt(stepRn)
                      << std::fixed
                      << "  loss=" << std::setprecision(6) << curLoss
                      << "  mean=" << std::setprecision(2) << curMeanMM << "mm"
                      << "  max=" << curMaxMM << "mm"
                      << (improved ? "  *" : "  (нет улучшения — стоп)") << "\n";
            if (!improved || evalNo >= budget) break;
        }
        std::cout << "\n========== ИТОГ (SVD / response-matrix) ==========\n"
                  << "  метод: Гаусс-Ньютон с усечённым SVD (коррекция орбиты)\n"
                  << "  бюджет: " << evalNo << " симуляций\n"
                  << std::fixed << std::setprecision(3)
                  << "  mean:  " << gBeforeMean << " mm  ->  " << curMeanMM << " mm\n"
                  << "  max:   " << gBeforeMax  << " mm  ->  " << curMaxMM  << " mm\n"
                  << "==================================================\n" << std::defaultfloat;
    } else if (lsGlobal) {
        // ====== GLOBAL LINE-SEARCH (последовательное приближение по всей конфигурации) ======
        // Покоординатный спуск: проходим по ВСЕМ группам, каждую настраиваем
        // детерминированным coarse-to-fine 1D-поиском по её допустимому диапазону
        // (loss глобальный, через runTrackingAndLoss). Зум вокруг лучшей точки
        // (включая incumbent), уполовинивание интервала. n раундов/элемент из бюджета
        // (cfg.nIter оценок) равномерно, остаток -> повторные свипы. dynamic_loss и
        // кики не применяются (метод детерминирован, монотонно не ухудшает).
        const int K   = std::max(3, cfg.perElement.lineSearchPoints);
        const int nG  = (int)rgroups.size();
        // ЕДИНЫЙ БЮДЖЕТ: total_iters (как в оконном режиме), fallback на n_iterations
        // для совместимости. Тот же распределитель budgetRounds, что везде.
        const int budget = std::max(1, cfg.perElement.totalIters > 0
                                          ? cfg.perElement.totalIters : cfg.nIter);
        const int blockSz = std::max(1, cfg.perElement.lineSearchBlockSize);
        const int nBlocks = (nG + blockSz - 1) / blockSz;
        // Раунды считаем на ЕДИНИЦУ line-search: по-группно это nG групп, блочно -- nBlocks
        // блоков (меньше единиц -> больше раундов зума на блок, иначе t только {-1,0,1} и
        // блок не двигается: ±полный радиус совместно расходит пучок).
        const int lsUnits = (blockSz <= 1) ? nG : nBlocks;
        const int roundsPerGroup = budgetRounds(budget, lsUnits, K, cfg.perElement.lineSearchMaxRounds, cfg.perElement.lineSearchMinRounds);
        // Стоимость свипа: по-группно nG·rounds·K; блочно -- nG проб градиента +
        // на каждый блок (incumbent + K·rounds line-search скаляра t вдоль -grad).
        const int sweepCost = (blockSz <= 1)
            ? std::max(1, nG * roundsPerGroup * K)
            : std::max(1, nG + nBlocks * (1 + K * roundsPerGroup));
        const int sweeps = std::max(1, budget / sweepCost);
        std::cout << "[GLOBAL LINE-SEARCH] groups=" << nG << "  block=" << blockSz
                  << "  K=" << K << "  rounds/grp=" << roundsPerGroup
                  << " (>=" << std::min(cfg.perElement.lineSearchMinRounds, cfg.perElement.lineSearchMaxRounds)
                  << ")  sweeps=" << sweeps
                  << "  диапазоны: |ток|<=" << cfg.perElement.currentRangeAmp << "A"
                  << ", |coeff|<=" << cfg.perElement.coeffRange
                  << "  (~" << sweeps * sweepCost << " оценок)\n";
        if (cfg.perElement.lineSearchLocalCurrent > 0.0 || cfg.perElement.lineSearchLocalCoeff > 0.0)
            std::cout << "  ЛОКАЛЬНЫЙ поиск: окрестность вокруг текущего значения  ток ±"
                      << cfg.perElement.lineSearchLocalCurrent << "A, coeff ±"
                      << cfg.perElement.lineSearchLocalCoeff << "  (против railing'а)\n";
        const double INF = std::numeric_limits<double>::infinity();
        int evalNo = 0;
        const double gBeforeMean = curMeanMM, gBeforeMax = curMaxMM;  // для ИТОГА
        // Печать ВСЕХ компонент loss (для подбора весов: при весах=1 видны сырые
        // величины членов, потом решаем какие веса дать).
        auto printComponents = [&](const std::string& label) {
            double L=0, mn=0, mx=0, cPos=0, cDir=0, cEnv=0, cPar=0, cDec=0;
            if (!runTrackingAndLoss(L, mn, mx, &cPos, &cDir, &cEnv, &cPar, &cDec)) {
                std::cout << "  [компоненты | " << label << "] пучок не дошёл до финиша\n";
                return;
            }
            std::cout << std::scientific << std::setprecision(4)
                      << "  [компоненты loss | " << label << "]  position=" << cPos
                      << "  parallelism=" << cPar
                      << "  envelope_decay=" << cDec << "\n"
                      << "      взвешенные:  pos=" << cPos
                      << "  par(w=" << std::defaultfloat << cfg.parallelismWeight << ")="
                      << std::scientific << cfg.parallelismWeight*cPar
                      << "  dec(w=" << std::defaultfloat << cfg.envelopeDecayWeight << ")="
                      << std::scientific << cfg.envelopeDecayWeight*cDec
                      << "  sum=" << L << std::defaultfloat << "\n";
        };
        printComponents("baseline");
        // Натуральный кубический сплайн через (x,y); возвращает x*, минимизирующий сплайн на
        // [x0, x_{n-1}] (плотная выборка 200 точек). Фоллбэк на argmin проб, если точек <3 или
        // среди y есть нефинитные (упавшая симуляция).
        auto splineMin = [](const std::vector<double>& x, const std::vector<double>& y,
                            double& yPred) -> double {
            const int n = (int)x.size();
            int jb = 0; for (int i = 1; i < n; ++i) if (y[i] < y[jb]) jb = i;
            bool fin = true; for (double v : y) if (!std::isfinite(v)) fin = false;
            if (n < 3 || !fin) { yPred = y[jb]; return x[jb]; }
            std::vector<double> aa(n,0.0), bb(n,0.0), cc(n,0.0), dd(n,0.0), mm(n,0.0);
            bb[0] = 1.0; bb[n-1] = 1.0;                      // natural: m[0]=m[n-1]=0
            for (int i = 1; i < n-1; ++i) {
                double h1 = x[i]-x[i-1], h2 = x[i+1]-x[i];
                aa[i] = h1; bb[i] = 2.0*(h1+h2); cc[i] = h2;
                dd[i] = 6.0*((y[i+1]-y[i])/h2 - (y[i]-y[i-1])/h1);
            }
            for (int i = 1; i < n; ++i) {                    // Thomas (прямой ход)
                double w = aa[i]/bb[i-1]; bb[i] -= w*cc[i-1]; dd[i] -= w*dd[i-1];
            }
            mm[n-1] = dd[n-1]/bb[n-1];
            for (int i = n-2; i >= 0; --i) mm[i] = (dd[i] - cc[i]*mm[i+1])/bb[i];
            const int S = 200; double best = std::numeric_limits<double>::infinity(), xb = x[jb];
            for (int s = 0; s <= S; ++s) {
                double xx = x[0] + (x[n-1]-x[0])*s/(double)S;
                int i = 0; while (i < n-2 && xx > x[i+1]) ++i;
                double h = x[i+1]-x[i], t = xx - x[i];
                double yy = y[i] + ((y[i+1]-y[i])/h - h*(2.0*mm[i]+mm[i+1])/6.0)*t
                                 + mm[i]*0.5*t*t + (mm[i+1]-mm[i])/(6.0*h)*t*t*t;
                if (yy < best) { best = yy; xb = xx; }
            }
            yPred = best;                                    // предсказание сплайна в минимуме
            return xb;
        };
        // МНК-парабола по ВСЕМ K точкам (least-squares y≈A·u²+B·u+C, u=x-mean — центрируем для
        // обусловленности). Усредняет дискретизац. шум (устойчивее 3-точечной). Вершина
        // x*=mean−B/(2A). Защита: A>0 (открыта вверх => минимум) и x* в диапазоне проб; иначе
        // откат на лучшую пробу. yPred = предсказание параболы в вершине (для лога).
        auto parabolaMin = [](const std::vector<double>& x, const std::vector<double>& y,
                              double& yPred) -> double {
            const int n = (int)x.size();
            int jb = 0; for (int i = 1; i < n; ++i) if (y[i] < y[jb]) jb = i;
            bool finAll = true; for (double v : y) if (!std::isfinite(v)) finAll = false;
            if (n < 3 || !finAll) { yPred = y[jb]; return x[jb]; }
            double xm = 0.0; for (double v : x) xm += v; xm /= n;
            // ВЕСА: loss по параметру — узкая яма, обычный МНК доминируется большими боковыми
            // значениями и мажет мимо. Взвешиваем в пользу низких loss => фит садится на дно ямы.
            double ymin = y[jb], ymax = y[0]; for (double v : y) ymax = std::max(ymax, v);
            const double delta = 0.1*(ymax - ymin) + 1e-12;
            double s0=0,s1=0,s2=0,s3=0,s4=0, t0=0,t1=0,t2=0;            // взвеш. Σw·uᵏ, Σw·uᵏy
            for (int i = 0; i < n; ++i) {
                double u = x[i]-xm, u2 = u*u;
                double w = 1.0/(y[i] - ymin + delta);
                s0+=w; s1+=w*u; s2+=w*u2; s3+=w*u2*u; s4+=w*u2*u2;
                t0+=w*y[i]; t1+=w*u*y[i]; t2+=w*u2*y[i];
            }
            // нормальные уравнения [s4 s3 s2; s3 s2 s1; s2 s1 s0]·[A;B;C]=[t2;t1;t0] (Cramer)
            const double D = s4*(s2*s0-s1*s1) - s3*(s3*s0-s1*s2) + s2*(s3*s1-s2*s2);
            if (std::abs(D) < 1e-300) { yPred = y[jb]; return x[jb]; }
            const double A = (t2*(s2*s0-s1*s1) - s3*(t1*s0-s1*t0) + s2*(t1*s1-s2*t0)) / D;
            const double B = (s4*(t1*s0-s1*t0) - t2*(s3*s0-s1*s2) + s2*(s3*t0-t1*s2)) / D;
            const double C = (s4*(s2*t0-t1*s1) - s3*(s3*t0-t1*s2) + t2*(s3*s1-s2*s2)) / D;
            if (A <= 1e-30) { yPred = y[jb]; return x[jb]; }            // не выпукла вверх -> нет min
            const double uStar = -B/(2.0*A);
            const double xStar = uStar + xm;
            double lo = x[0], hi = x[n-1]; if (lo > hi) std::swap(lo, hi);
            if (!(xStar > lo && xStar < hi)) { yPred = y[jb]; return x[jb]; }  // вне диапазона -> проба
            yPred = A*uStar*uStar + B*uStar + C;
            return xStar;
        };
        auto lsGroup = [&](int gi, int sweep) {
            auto& rg = rgroups[gi];
            const bool isCoeff = (rg.cfg.target.rfind("COEFF", 0) == 0);
            const double absLo = isCoeff ? -cfg.perElement.coeffRange : -cfg.perElement.currentRangeAmp;
            const double absHi = isCoeff ?  cfg.perElement.coeffRange :  cfg.perElement.currentRangeAmp;
            // Локальный радиус: если задан, ищем в ОКРЕСТНОСТИ текущего значения
            // [cur-R, cur+R] (зажатой в [absLo,absHi]) -- против railing'а.
            const double localR = isCoeff ? cfg.perElement.lineSearchLocalCoeff
                                          : cfg.perElement.lineSearchLocalCurrent;
            const double lo = (localR > 0.0) ? std::max(absLo, rg.currentValue - localR) : absLo;
            const double hi = (localR > 0.0) ? std::min(absHi, rg.currentValue + localR) : absHi;
            double bestVal = std::min(std::max(rg.currentValue, lo), hi);
            // incumbent (текущее значение элемента) -- среди кандидатов => монотонность
            applyGroupValue(rg, bestVal);
            double L = 0, mean = 0, mx = 0;
            bool ok = runTrackingAndLoss(L, mean, mx);
            ++evalNo;
            double bestLossG = ok ? L : INF;
            double bMean = ok ? mean : curMeanMM, bMax = ok ? mx : curMaxMM;
            if (ok && L < bestLoss) { bestLoss = L;
                for (int k = 0; k < nG; ++k) bestValues[k] = rgroups[k].currentValue; }
            const double incumbentLoss = bestLossG;
            auto tryX = [&](double x) -> double {        // оценить loss в x, обновить лучшее, вернуть loss
                applyGroupValue(rg, x);
                double Lx = 0, mnx = 0, mxx = 0;
                const bool okx = runTrackingAndLoss(Lx, mnx, mxx);
                ++evalNo;
                if (!okx) return INF;
                if (Lx < bestLossG) { bestLossG = Lx; bestVal = x; bMean = mnx; bMax = mxx; }
                if (Lx < bestLoss) { bestLoss = Lx;
                    for (int kk = 0; kk < nG; ++kk) bestValues[kk] = rgroups[kk].currentValue; }
                return Lx;
            };
            if (cfg.perElement.lineSearchParabola || cfg.perElement.lineSearchSpline) {
                // K(>=5) проб по [lo,hi]; аппроксимация loss(value) и её минимум на отрезке как
                // кандидат. PARABOLA (Брент) — парабола по тройке вокруг лучшей; SPLINE — кубич.
                // сплайн. Кандидат ВСЕГДА проверяется симуляцией => принимаем лучшую ОЦЕНЁННУЮ точку.
                const bool useParab = cfg.perElement.lineSearchParabola;   // приоритет у параболы
                // Полный ОДЗ (localR<=0) => диапазон широкий, нужно больше проб, иначе узкая яма
                // провалится между узлами. Локально (±R) хватает 5.
                const int Ksp = (localR > 0.0) ? std::max(5, K) : std::max(11, K);
                std::vector<double> xs(Ksp), ys(Ksp);
                for (int k = 0; k < Ksp; ++k) {
                    const double x = lo + (hi - lo) * k / (double)(Ksp - 1);
                    xs[k] = x; ys[k] = tryX(x);
                }
                double yPred = 0.0;
                const double xstar = useParab ? parabolaMin(xs, ys, yPred) : splineMin(xs, ys, yPred);
                const double Lact = tryX(xstar);         // кандидат -> проверить симуляцией
                // ЛОГ точек аппроксимации (для анализа модели). yhat = предсказание модели в x*;
                // Lact = реальный loss там. yhat<<Lact => модель переоценила провал (шум/осцилляция).
                std::cout << "  [" << (useParab ? "PARAB" : "SPLINE") << " g" << (gi+1)
                          << " " << rg.cfg.name << "] x=" << std::fixed << std::setprecision(4);
                for (int k = 0; k < Ksp; ++k) std::cout << (k ? "," : "") << xs[k];
                std::cout << " loss=" << std::scientific << std::setprecision(5);
                for (int k = 0; k < Ksp; ++k) std::cout << (k ? "," : "") << ys[k];
                std::cout << "  x*=" << std::fixed << std::setprecision(4) << xstar
                          << " yhat=" << std::scientific << std::setprecision(5) << yPred
                          << " Lact=" << Lact << std::defaultfloat << "\n";
            } else {
                double a = lo, b = hi;
                bool fiBreak = false;                    // first-improvement: ранний выход
                for (int r = 0; r < roundsPerGroup && !fiBreak; ++r) {
                    const double h = (b - a) / (K - 1);
                    for (int k = 0; k < K; ++k) {
                        const double x = a + k * h;
                        applyGroupValue(rg, x);
                        double Lx = 0, mnx = 0, mxx = 0;
                        const bool okx = runTrackingAndLoss(Lx, mnx, mxx);
                        ++evalNo;
                        if (okx && Lx < bestLossG) { bestLossG = Lx; bestVal = x; bMean = mnx; bMax = mxx; }
                        if (okx && Lx < bestLoss) { bestLoss = Lx;
                            for (int kk = 0; kk < nG; ++kk) bestValues[kk] = rgroups[kk].currentValue; }
                        // FIRST-IMPROVEMENT: первое значение, улучшившее тюн группы -> принять и
                        // НЕ дозумивать (экономит сим., грубее => легче выходит из острых экстремумов).
                        if (cfg.perElement.lineSearchFirstImprove && okx && Lx < incumbentLoss) {
                            fiBreak = true; break;
                        }
                    }
                    if (fiBreak) break;
                    const double newHalf = h * 0.5;      // уполовинивание интервала
                    a = std::max(lo, bestVal - newHalf);
                    b = std::min(hi, bestVal + newHalf);
                }
            }
            applyGroupValue(rg, bestVal);                // зафиксировать лучшее значение
            curLoss = bestLossG; curMeanMM = bMean; curMaxMM = bMax;
            const bool improved = (bestLossG < incumbentLoss);
            if (improved) accepted++;
            history.push_back(HistEntry{evalNo, rg.cfg.name, 0.0, bestVal, bestLossG,
                                        bestLoss, bMean, bMax, improved ? "LS+IMPROVE" : "LS"});
            std::cout << "[LS s" << sweep << "/" << sweeps
                      << " g" << (gi + 1) << "/" << nG << " "
                      << std::left << std::setw(24) << rg.cfg.name << std::right
                      << "] val=" << std::fixed << std::setprecision(4) << bestVal
                      << "  loss=" << std::setprecision(6) << bestLossG
                      << "  best=" << bestLoss
                      << "  mean=" << std::setprecision(2) << bMean << "mm"
                      << "  max=" << bMax << "mm"
                      << (improved ? "  *" : "") << "\n";
            if (evalNo % CHECKPOINT_EVERY == 0) writeHistory(false);
        };
        // БЛОЧНЫЙ ГРАДИЕНТНЫЙ line-search: blockSize групп оптимизируются ВМЕСТЕ.
        // На блок: градиент в u-пространстве (u=(x-x0)/локрадиус, конечн. разность du=0.1)
        // + line-search скаляра t∈[-1,1] вдоль -grad (нормировано => в пределах локрадиуса).
        // Двигает СВЯЗАННЫЕ токи совместно -> выходит из координатных застреваний.
        auto lsBlock = [&](int gStart, int gCount, int sweep) {
            std::vector<double> x0(gCount), locR(gCount), aLo(gCount), aHi(gCount), dirU(gCount, 0.0);
            for (int b = 0; b < gCount; ++b) {
                auto& rg = rgroups[gStart + b];
                const bool isC = (rg.cfg.target.rfind("COEFF", 0) == 0);
                double R = isC ? cfg.perElement.lineSearchLocalCoeff : cfg.perElement.lineSearchLocalCurrent;
                if (R <= 0.0) R = isC ? cfg.perElement.coeffRange : cfg.perElement.currentRangeAmp;
                locR[b] = R;
                aLo[b]  = isC ? -cfg.perElement.coeffRange : -cfg.perElement.currentRangeAmp;
                aHi[b]  = isC ?  cfg.perElement.coeffRange :  cfg.perElement.currentRangeAmp;
                x0[b]   = rg.currentValue;
            }
            auto setT = [&](double t) {
                for (int b = 0; b < gCount; ++b) {
                    double x = x0[b] + t * dirU[b] * locR[b];
                    applyGroupValue(rgroups[gStart + b], std::min(aHi[b], std::max(aLo[b], x)));
                }
            };
            double L0 = 0, m0 = 0, mx0 = 0;                       // incumbent (t=0)
            bool ok0 = runTrackingAndLoss(L0, m0, mx0); ++evalNo;
            if (!ok0) L0 = INF;
            if (ok0 && L0 < bestLoss) { bestLoss = L0;
                for (int k = 0; k < nG; ++k) bestValues[k] = rgroups[k].currentValue; }
            // du мал (0.02 локрадиуса) -> проба остаётся в ГЛАДКОЙ зоне, не прыгает за
            // обрыв расходимости. Кап Lc на удвоение L0 не даёт «обрыву» (Lp огромный
            // или окп=false) доминировать в направлении -- худшая ось даёт ограниченный
            // штраф, а не 1e6. Улучшающие пробы (Lp<L0) сохраняют полный градиент.
            const double du = 0.02;
            const double Lcap = (L0 < INF ? L0 : 1.0) * 2.0 + 1e-9;
            std::vector<double> gU(gCount, 0.0);
            int nFail = 0;
            for (int b = 0; b < gCount; ++b) {
                auto& rg = rgroups[gStart + b];
                double xp = std::min(aHi[b], std::max(aLo[b], x0[b] + du * locR[b]));
                applyGroupValue(rg, xp);
                double Lp = 0, mp = 0, mxp = 0;
                bool okp = runTrackingAndLoss(Lp, mp, mxp); ++evalNo;
                applyGroupValue(rg, x0[b]);
                double aDu = (xp - x0[b]) / locR[b];
                double Lc = okp ? std::min(Lp, Lcap) : Lcap;
                gU[b] = (std::abs(aDu) > 1e-9) ? (Lc - L0) / aDu : 0.0;
                if (!okp) ++nFail;
            }
            double gn = 0.0; for (int b = 0; b < gCount; ++b) gn += gU[b]*gU[b];
            gn = std::sqrt(gn);
            for (int b = 0; b < gCount; ++b) dirU[b] = (gn > 1e-30) ? -gU[b] / gn : 0.0;
            double bestT = 0.0, bestLossB = L0, bMean = m0, bMax = mx0;
            if (gn > 1e-30) {                                    // line-search t вдоль -grad
                double aT = -1.0, bT = 1.0;
                for (int r = 0; r < roundsPerGroup; ++r) {
                    const double h = (bT - aT) / (K - 1);
                    for (int k = 0; k < K; ++k) {
                        const double t = aT + k * h;
                        setT(t);
                        double Lt = 0, mt = 0, mxt = 0;
                        bool okt = runTrackingAndLoss(Lt, mt, mxt); ++evalNo;
                        if (okt && Lt < bestLossB) { bestLossB = Lt; bestT = t; bMean = mt; bMax = mxt; }
                        if (okt && Lt < bestLoss) { bestLoss = Lt;
                            for (int kk = 0; kk < nG; ++kk) bestValues[kk] = rgroups[kk].currentValue; }
                    }
                    const double nh = h * 0.5;
                    aT = std::max(-1.0, bestT - nh);
                    bT = std::min( 1.0, bestT + nh);
                }
            }
            setT(bestT);                                         // зафиксировать лучшее (t=0 => x0)
            curLoss = bestLossB; curMeanMM = bMean; curMaxMM = bMax;
            const bool improved = (bestLossB < L0);
            if (improved) accepted++;
            history.push_back(HistEntry{evalNo, "BLOCK", 0.0, bestT, bestLossB,
                                        bestLoss, bMean, bMax, improved ? "LSB+IMPROVE" : "LSB"});
            std::cout << "[LSB s" << sweep << "/" << sweeps << " блок "
                      << (gStart + 1) << "-" << (gStart + gCount) << "/" << nG
                      << "] t=" << std::fixed << std::setprecision(3) << bestT
                      << "  loss=" << std::setprecision(6) << bestLossB
                      << "  best=" << bestLoss
                      << "  grad=" << std::scientific << std::setprecision(2) << gn
                      << " fail=" << nFail << "/" << gCount << std::fixed
                      << "  max=" << std::setprecision(2) << bMax << "mm"
                      << (improved ? "  *" : "") << "\n";
            if (evalNo % CHECKPOINT_EVERY == 0) writeHistory(false);
        };
        // first-improve удешевляет свип => даём ему делать БОЛЬШЕ свипов до исчерпания бюджета
        // (грубые повторные проходы лучше борются с экстремумами, чем один точный). Иначе
        // фиксированное число свипов как раньше.
        // first-improve и сплайн дешевле полного зума => даём им крутить БОЛЬШЕ свипов до
        // исчерпания бюджета (стоп по бюджету, в т.ч. мид-свип). Иначе фикс. число свипов.
        const bool budgetBounded = cfg.perElement.lineSearchFirstImprove
                                 || cfg.perElement.lineSearchSpline || cfg.perElement.lineSearchParabola;
        const int maxSweeps = budgetBounded ? std::max(sweeps, budget) : sweeps;
        const bool fiBudget = budgetBounded;
        for (int s = 1; s <= maxSweeps; ++s) {
            if (blockSz <= 1)
                for (int gi = 0; gi < nG; ++gi) {
                    lsGroup(gi, s);
                    if (fiBudget && evalNo >= budget) break;          // мид-свип стоп => без перерасхода
                }
            else
                for (int gStart = 0; gStart < nG; gStart += blockSz) {
                    lsBlock(gStart, std::min(blockSz, nG - gStart), s);
                    if (fiBudget && evalNo >= budget) break;
                }
            // Компоненты ЛУЧШЕГО состояния после свипа (без изменения хода поиска:
            // временно применяем bestValues, печатаем, возвращаем текущие).
            std::vector<double> cur(nG);
            for (int k = 0; k < nG; ++k) { cur[k] = rgroups[k].currentValue; applyGroupValue(rgroups[k], bestValues[k]); }
            printComponents("sweep " + std::to_string(s));
            for (int k = 0; k < nG; ++k) applyGroupValue(rgroups[k], cur[k]);
            // КИК: каждые kickEverySweeps свипов встряхиваем ВСЕ группы на ±amp×|value|
            // (только если впереди есть свип для доводки). bestValues не трогаем — плохой
            // кик не потеряет лучшее; следующий свип доводит из встряхнутого состояния.
            if (cfg.perElement.lineSearchKickEverySweeps > 0 && s < maxSweeps
                    && s % cfg.perElement.lineSearchKickEverySweeps == 0) {
                const double amp = cfg.perElement.lineSearchKickAmplitude;
                for (auto& rg : rgroups) {
                    const bool isC = (rg.cfg.target.rfind("COEFF", 0) == 0);
                    const double absHi = isC ? cfg.perElement.coeffRange : cfg.perElement.currentRangeAmp;
                    double nv = rg.currentValue + rg.currentValue * amp * uniMP(rng);
                    applyGroupValue(rg, std::min(std::max(nv, -absHi), absHi));
                }
                double km = 0, kx = 0, kl = 0;
                if (runTrackingAndLoss(kl, km, kx)) { curLoss = kl; curMeanMM = km; curMaxMM = kx; }
                ++evalNo;
                std::cout << "[LS КИК после свипа " << s << "]  ±"
                          << std::fixed << std::setprecision(1) << (amp * 100.0)
                          << "% от значения  loss=" << std::setprecision(6) << curLoss
                          << "  best=" << bestLoss << "\n";
            }
            if (fiBudget && evalNo >= budget) break;                  // бюджет исчерпан
        }
        // ============== ИТОГ (глобальная оптимизация) ==============
        {
            auto pct = [](double a, double b) {
                return (std::abs(a) > 1e-30) ? 100.0 * (a - b) / a : 0.0; };
            std::cout << "\n========== ИТОГ (глобальная оптимизация) ==========\n";
            std::cout << "  метод: line-search (последовательное приближение по всей конфигурации)";
            if (cfg.perElement.lineSearchLocalCurrent > 0.0)
                std::cout << ", локальный радиус ±" << std::fixed << std::setprecision(0)
                          << cfg.perElement.lineSearchLocalCurrent << "A";
            std::cout << "\n  бюджет: " << evalNo << " симуляций, " << sweeps << " свипов\n";
            std::cout << std::fixed << std::setprecision(3)
                      << "  mean:       " << gBeforeMean << " mm  ->  " << curMeanMM
                      << " mm    (" << std::setprecision(2) << pct(gBeforeMean, curMeanMM) << "%)\n"
                      << std::setprecision(3)
                      << "  max:        " << gBeforeMax << " mm  ->  " << curMaxMM
                      << " mm    (" << std::setprecision(2) << pct(gBeforeMax, curMaxMM) << "%)\n";
            std::vector<double> absI;
            for (const auto& rg : rgroups)
                if (rg.cfg.target.rfind("COEFF", 0) != 0 && std::abs(rg.currentValue) > 1e-9)
                    absI.push_back(std::abs(rg.currentValue));
            if (!absI.empty()) {
                std::sort(absI.begin(), absI.end());
                int railed = 0;
                for (double v : absI) if (v >= cfg.perElement.currentRangeAmp - 0.1) ++railed;
                std::cout << std::fixed << std::setprecision(1)
                          << "  токи:       " << absI.size() << " ненул., медиана|I|="
                          << absI[absI.size()/2] << "A, max|I|=" << absI.back()
                          << "A, на упоре=" << railed << "\n";
            }
            std::cout << "  время:      " << std::fixed << std::setprecision(1)
                      << totalTimer.sec() << "s\n" << std::defaultfloat;
            std::cout << "==================================================\n";
        }
    } else {
    std::cout << "[ГЛОБАЛЬНЫЙ РЕЖИМ] метод: monte_carlo (случайные кики, hill-climber)"
              << "  n_iter=" << cfg.nIter;
    if (kickTrigger > 0) std::cout << "  кики каждые " << kickTrigger;
    std::cout << "\n";
    for (int it = 1; it <= cfg.nIter; ++it) {
        bool isKick = (kickTrigger > 0 && itersSinceKick >= kickTrigger);

        if (isKick) {
            // KICK: возмущаем все, безусловно принимаем
            std::vector<double> oldVals(rgroups.size());
            for (size_t k = 0; k < rgroups.size(); ++k) {
                oldVals[k] = rgroups[k].currentValue;
                double d = randomStep(rgroups[k], cfg.kickScale);
                applyGroupValue(rgroups[k], oldVals[k] + d);
            }
            double newLoss = 0, meanMM = 0, maxMM = 0;
            if (!runTrackingAndLoss(newLoss, meanMM, maxMM)) {
                // Откат
                for (size_t k = 0; k < rgroups.size(); ++k) applyGroupValue(rgroups[k], oldVals[k]);
                failed++;
                HistEntry e{it, "KICK", 0, 0, curLoss, bestLoss, curMeanMM, curMaxMM, "KICK_FAIL"};
                history.push_back(e);
                std::cout << "[" << it << "/" << cfg.nIter << "] KICK FAILED, rolled back"
                          << "  best=" << bestLoss << "\n";
                continue;
            }
            curLoss = newLoss;
            curMeanMM = meanMM; curMaxMM = maxMM;
            std::string mark = "KICK";
            if (newLoss < bestLoss) {
                bestLoss = newLoss;
                for (size_t k = 0; k < rgroups.size(); ++k) bestValues[k] = rgroups[k].currentValue;
                mark = "KICK+BEST";
                kicksWithoutBest = 0;
            } else {
                kicksWithoutBest++;
            }
            kicks++;
            itersSinceKick = 0;
            HistEntry e{it, "KICK", 0, 0, newLoss, bestLoss, meanMM, maxMM, mark};
            history.push_back(e);
            std::cout << "[" << it << "/" << cfg.nIter << "] *** " << mark
                      << "  loss=" << std::fixed << std::setprecision(6) << newLoss
                      << "  best=" << bestLoss
                      << "  mean=" << std::setprecision(2) << meanMM << "mm"
                      << "  max=" << maxMM << "mm"
                      << " (kicks-no-best=" << kicksWithoutBest << ")\n";

            // Если N кикcов подряд не дали best -- откат к best как стартовой точке.
            // При REVERT_AFTER_KICKS=0 функция отключена.
            if (REVERT_AFTER_KICKS > 0 && kicksWithoutBest >= REVERT_AFTER_KICKS) {
                for (size_t k = 0; k < rgroups.size(); ++k)
                    applyGroupValue(rgroups[k], bestValues[k]);
                // Пересчитаем loss с best-параметрами для cur_loss
                double mm = 0, mx = 0, newL = 0;
                if (runTrackingAndLoss(newL, mm, mx)) {
                    curLoss = newL; curMeanMM = mm; curMaxMM = mx;
                }
                reverts++;
                kicksWithoutBest = 0;
                itersSinceKick = 0;   // сбросим, чтобы кик не сработал сразу же
                HistEntry er{it, "REVERT", 0, 0, curLoss, bestLoss, curMeanMM, curMaxMM, "REVERT"};
                history.push_back(er);
                std::cout << "[" << it << "/" << cfg.nIter << "] *** REVERT to best"
                          << "  cur_loss=" << std::fixed << std::setprecision(6) << curLoss
                          << "  best=" << bestLoss << "\n";
            }
            continue;
        }

        // Normal iteration: одна группа, random step, greedy accept
        RuntimeGroup& rg = rgroups[normalIdx % rgroups.size()];
        normalIdx++;
        itersSinceKick++;
        double oldVal = rg.currentValue;
        double delta = randomStep(rg, 1.0);
        double newVal = oldVal + delta;
        applyGroupValue(rg, newVal);

        double newLoss = 0, meanMM = 0, maxMM = 0;
        // Запрашиваем компоненты unweighted -- они нужны и для dynamic_loss, и
        // для записи в history. При выключенном dynamic_loss оверхед ~5%, не критично.
        double new_pos_c = 0, new_dir_c = 0, new_env_c = 0;
        if (!runTrackingAndLoss(newLoss, meanMM, maxMM,
                                 &new_pos_c, &new_dir_c, &new_env_c)) {
            applyGroupValue(rg, oldVal);
            failed++;
            HistEntry e{it, rg.cfg.name, oldVal, newVal, curLoss, bestLoss, curMeanMM, curMaxMM, "FAIL"};
            history.push_back(e);
            std::cout << "[" << it << "/" << cfg.nIter << "] " << rg.cfg.name
                      << "  " << oldVal << " -> " << newVal << "  FAIL  best=" << bestLoss << "\n";
            continue;
        }
        // Решение о accept: при dynamic_loss enabled сравниваем по L_dyn(t),
        // иначе -- по статическому composite (legacy).
        bool accept_step;
        if (cfg.dynamicLoss.enabled) {
            double wp, wd, we;
            dynamicWeights(it, wp, wd, we);
            double L_cur_dyn = composeLoss(wp, wd, we, cur_pos_c, cur_dir_c, cur_env_c);
            double L_new_dyn = composeLoss(wp, wd, we, new_pos_c, new_dir_c, new_env_c);
            accept_step = (L_new_dyn < L_cur_dyn);
        } else {
            accept_step = (newLoss < curLoss);
        }
        std::string mark;
        if (accept_step) {
            curLoss = newLoss;
            curMeanMM = meanMM; curMaxMM = maxMM;
            // Обновляем компоненты cur_*: при следующей итерации они будут
            // использоваться как baseline для динамического сравнения.
            cur_pos_c = new_pos_c; cur_dir_c = new_dir_c; cur_env_c = new_env_c;
            accepted++;
            mark = "ACCEPT";
            if (newLoss < bestLoss) {
                bestLoss = newLoss;
                for (size_t k = 0; k < rgroups.size(); ++k) bestValues[k] = rgroups[k].currentValue;
                mark = "ACCEPT+BEST";
                kicksWithoutBest = 0;  // сбросим счётчик при любом новом best
            }
        } else {
            applyGroupValue(rg, oldVal);
            rejected++;
            mark = "reject";
        }
        HistEntry e{it, rg.cfg.name, oldVal, newVal, newLoss, bestLoss, meanMM, maxMM, mark};
        history.push_back(e);

        // Логи -- редкие, чтобы не тормозить вывод
        if (it <= 20 || it % 50 == 0 || mark.find("BEST") != std::string::npos) {
            std::cout << "[" << it << "/" << cfg.nIter << "] " << std::left << std::setw(28) << rg.cfg.name
                      << " " << std::right << std::setw(12) << std::fixed << std::setprecision(6) << oldVal
                      << " -> " << std::setw(12) << newVal
                      << "  loss=" << std::setprecision(6) << newLoss
                      << "  best=" << bestLoss
                      << "  (mean=" << std::setprecision(2) << meanMM << "mm"
                      << " max=" << maxMM << "mm)"
                      << "  [" << mark << "]\n";
        }

        // Периодический чекпойнт истории на случай Ctrl+C / краша.
        if (it % CHECKPOINT_EVERY == 0) {
            writeHistory(false);
        }
    }
    }  // end else (random global mode)
    }  // end of MC loop wrapper {} (dynamic_loss scope)
skip_mc_loop:
    double loopSec = loopTimer.sec();

    // ---- 13. Финализация: BEST или (при accept-final) ТЕКУЩЕЕ исследованное ----
    // accept-final разрешает «откат» (сохранить состояние хуже best, чтобы из него
    // шагнуть дальше) -- но НЕ катастрофический разлёт: если loss >> best, безопасно
    // откатываемся к best.
    {
        const double ROLLBACK_CAP = 50.0;
        const bool keepFinal = cfg.perElement.lineSearchAcceptFinal
                               && std::isfinite(curLoss) && curLoss < bestLoss * ROLLBACK_CAP;
        if (keepFinal) {
            std::cout << "[accept-final] оставлено ТЕКУЩЕЕ состояние (loss="
                      << std::scientific << std::setprecision(4) << curLoss
                      << ", best=" << bestLoss << std::defaultfloat << ") -- откат разрешён\n";
        } else {
            for (size_t k = 0; k < rgroups.size(); ++k) applyGroupValue(rgroups[k], bestValues[k]);
        }
    }

    // ---- 13b. Финальный прогон + запись trajectories_exp{id}.bin ----
    // Иначе visualize.py покажет СТАРУЮ траекторию из предыдущего запуска
    // SynchrotronTracker.exe -- что вводит в заблуждение.
    {
        std::vector<std::vector<FieldVector>> finalTraj;
        if (runTracking(finalTraj)) {
            std::string binP = ws + "results/trajectories_exp" + std::to_string(cfg.expId) + ".bin";
            std::ofstream out(binP, std::ios::binary);
            for (size_t i = 0; i < finalTraj.size(); ++i) {
                uint32_t sz = (uint32_t)finalTraj[i].size();
                out.write((char*)&sz, 4);
                out.write((char*)finalTraj[i].data(), sz * sizeof(FieldVector));
            }
            out.close();
            std::cout << "Final trajectory: " << binP
                      << " (" << finalTraj.size() << " particles)\n";
        } else {
            std::cerr << "WARN: final tracking failed; trajectory bin not updated\n";
        }
    }

    // ---- 14. Save history JSON (финальная запись) ----
    writeHistory(true);
    std::cout << "\nHistory: " << histPath << "\n";

    // ---- 15. Save best state back to CSVs ----
    int sp_idx_val_w   = spCsv.col("arg_val");
    int sp_idx_coeff_w = spCsv.col("correct_coeff");
    bool spDirty = false;
    for (auto& rg : rgroups) {
        if (rg.cfg.target == "CURRENT") {
            for (auto& key : rg.spKeys) {
                auto it = spMap.find(key);
                if (it == spMap.end()) continue;
                auto& row = spCsv.rows[it->second.row_index];
                if ((int)row.size() > sp_idx_val_w) {
                    std::ostringstream os;
                    os << std::fixed << std::setprecision(10) << rg.currentValue;
                    row[sp_idx_val_w] = os.str();
                    spDirty = true;
                }
            }
        } else if (rg.cfg.target == "COEFF_INST" && sp_idx_coeff_w >= 0) {
            for (auto& key : rg.spKeys) {
                auto it = spMap.find(key);
                if (it == spMap.end()) continue;
                auto& row = spCsv.rows[it->second.row_index];
                if ((int)row.size() > sp_idx_coeff_w) {
                    std::ostringstream os;
                    os << std::fixed << std::setprecision(10) << rg.currentValue;
                    row[sp_idx_coeff_w] = os.str();
                    spDirty = true;
                }
            }
        }
    }
    if (spDirty) {
        writeCsvTable(ws + "tables/sp_elements.csv", spCsv);
        std::cout << "Saved best -> sp_elements.csv\n";
    }
    // Также обновим hardware_library для COEFF_HW-групп
    int hw_cc_w = hwCsv.col("correct_coeff");
    bool hwDirty = false;
    for (auto& rg : rgroups) {
        if (rg.cfg.target == "COEFF_HW" && hw_cc_w >= 0) {
            for (auto& dt : rg.hwDeviceTypes) {
                for (auto& r : hwCsv.rows) {
                    if ((int)r.size() > hw_id && r[hw_id] == dt && (int)r.size() > hw_cc_w) {
                        std::ostringstream os;
                        os << std::fixed << std::setprecision(10) << rg.currentValue;
                        r[hw_cc_w] = os.str();
                        hwDirty = true;
                    }
                }
            }
        }
    }
    if (hwDirty) {
        writeCsvTable(ws + "tables/hardware_library.csv", hwCsv);
        std::cout << "Saved best -> hardware_library.csv\n";
    }

    // Пересборка lattice_configs.csv для совместимости с SynchrotronTracker.exe.
    // ВАЖНО: при пересборке мы пишем ТОЛЬКО строки, относящиеся к текущей
    // lattice_id (latId). Если в файле есть другие lattice_id -- их сохраняем
    // как есть (читаем оригинал перед перезаписью). И применяем latticeMaxGid.
    {
        // 1. Считать существующий lattice_configs.csv и сохранить чужие lattice_id.
        std::vector<std::vector<std::string>> preservedRows;
        {
            auto existing = readCsvTable(ws + "tables/lattice_configs.csv");
            for (const auto& r : existing.rows) {
                if (r.empty()) continue;
                int rowLat;
                try { rowLat = std::stoi(r[0]); } catch (...) { continue; }
                if (rowLat != latId) preservedRows.push_back(r);
            }
        }

        CsvTable lc;
        lc.header = {"lattice_id", "instance_id", "device_type_id",
                     "x", "y", "z", "yaw", "pitch", "roll",
                     "arg_val", "correct_coeff"};
        // Сначала -- строки чужих lattice (сохраняем их).
        for (const auto& r : preservedRows) lc.rows.push_back(r);
        for (auto& r : layoutCsv.rows) {
            if ((int)r.size() <= lay_roll) continue;
            // Применяем тот же фильтр по global_id, что и при сборке lattice.
            if (cfg.latticeMaxGid > 0) {
                int gid;
                try { gid = std::stoi(r[lay_gid]); } catch (...) { continue; }
                if (gid > cfg.latticeMaxGid) continue;
            }
            std::string spt = r[lay_typ];
            int seid;
            try { seid = std::stoi(r[lay_eid]); } catch (...) { continue; }
            auto it = spMap.find({spt, seid});
            if (it == spMap.end()) continue;
            // Берём актуальный arg_val и correct_coeff после правок.
            const auto& srow = spCsv.rows[it->second.row_index];
            std::string argStr   = srow[sp_idx_val_w];
            std::string coeffStr = (sp_idx_coeff_w >= 0 && (int)srow.size() > sp_idx_coeff_w)
                                   ? srow[sp_idx_coeff_w] : "1.0";
            std::vector<std::string> outRow = {
                std::to_string(latId),
                r[lay_gid],
                it->second.device_type_id,
                r[lay_x], r[lay_y], r[lay_z],
                r[lay_yaw], r[lay_pit], r[lay_roll],
                argStr, coeffStr,
            };
            lc.rows.push_back(std::move(outRow));
        }
        writeCsvTable(ws + "tables/lattice_configs.csv", lc);
        std::cout << "Rebuilt lattice_configs.csv (" << lc.rows.size() << " rows)\n";
    }

    std::cout << "\nDone in " << std::fixed << std::setprecision(1) << loopSec << "s: "
              << accepted << " accepted, " << rejected << " rejected, "
              << kicks << " kicks, " << reverts << " reverts, "
              << failed << " failed\n";
    std::cout << "Best loss: " << std::setprecision(6) << bestLoss << "\n";
    std::cout << "Total: " << totalTimer.sec() << "s\n";
    return 0;
}
