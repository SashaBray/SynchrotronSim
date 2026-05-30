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
    std::vector<int> lossSpNumbers;
    // Если > 0, в lattice включаем только элементы lattice_layout с
    // global_id <= latticeMaxGid (для экспериментов с "обрезанной" lattice,
    // например пролёт через первые 3 SP с последующим свободным дрейфом).
    int latticeMaxGid = 0;
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
    std::string historyOut = "optim_history.json";
    std::vector<GroupConfig> groups;
};

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
            if      (key == "EXP_ID")        cfg.expId      = std::stoi(val);
            else if (key == "N_ITER")        cfg.nIter      = std::stoi(val);
            else if (key == "SEED")          cfg.seed       = std::stoi(val);
            else if (key == "KICK_EVERY")    cfg.kickEvery  = std::stoi(val);
            else if (key == "KICK_SCALE")    cfg.kickScale  = std::stod(val);
            else if (key == "REVERT_AFTER_KICKS") cfg.revertAfterKicks = std::stoi(val);
            else if (key == "LOSS_TYPE")     cfg.lossType   = val;
            else if (key == "LOSS_SQUARED")  cfg.lossSquared= (std::stoi(val) != 0);
            else if (key == "LOSS_SP") {
                cfg.lossSpNumbers.clear();
                std::stringstream ss(val);
                int v;
                while (ss >> v) cfg.lossSpNumbers.push_back(v);
            }
            else if (key == "LATTICE_MAX_GID") cfg.latticeMaxGid = std::stoi(val);
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
            else if (key == "LOSS_EXCLUDE_DEVICES") {
                cfg.lossExcludeDevices.clear();
                std::stringstream ss(val);
                std::string dev;
                while (ss >> dev) cfg.lossExcludeDevices.push_back(dev);
            }
            else if (key == "HISTORY_OUT")   cfg.historyOut = val;
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
        std::cout << "Usage: SynchrotronOptimizer.exe <config.optimconfig>\n";
        std::cout << "  Сначала запустите prepare_optim_config.py для конвертации JSON -> optimconfig.\n";
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
    std::cout << "Config: " << cfgPath << "\n";
    std::cout << "  exp_id=" << cfg.expId
              << "  n_iter=" << cfg.nIter
              << "  kick_every=" << cfg.kickEvery
              << "  kick_scale=" << cfg.kickScale
              << "  groups=" << cfg.groups.size() << "\n";

    // ---- 2. Load experiment ----
    auto expCsv = readCsvTable(ws + "tables/experiments.csv");
    int latId = -1, grpId = -1, maxTurns = 0;
    int trajStride = 1;
    double dt = 0.0;
    FieldVector fp1 = {5.0, -1.0, -1.0};
    FieldVector fp2 = {5.0,  1.0, -1.0};
    FieldVector fp3 = {5.0,  0.0,  1.0};
    for (auto& r : expCsv.rows) {
        if ((int)r.size() < 5) continue;
        if (std::stoi(r[0]) == cfg.expId) {
            latId    = std::stoi(r[1]);
            grpId    = std::stoi(r[2]);
            dt       = std::stod(r[3]);
            maxTurns = std::stoi(r[4]);
            // Колонка 5 -- traj_stride. Backward compat: если в строке 14
            // столбцов -- старый формат без traj_stride.
            if (r.size() >= 15) {
                try { trajStride = std::max(1, std::stoi(r[5])); } catch (...) {}
                fp1 = {std::stod(r[6]),  std::stod(r[7]),  std::stod(r[8])};
                fp2 = {std::stod(r[9]),  std::stod(r[10]), std::stod(r[11])};
                fp3 = {std::stod(r[12]), std::stod(r[13]), std::stod(r[14])};
            } else if (r.size() >= 14) {
                fp1 = {std::stod(r[5]),  std::stod(r[6]),  std::stod(r[7])};
                fp2 = {std::stod(r[8]),  std::stod(r[9]),  std::stod(r[10])};
                fp3 = {std::stod(r[11]), std::stod(r[12]), std::stod(r[13])};
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
            e.arg_val        = (r[sp_idx_val] == "-") ? 0.0 : std::stod(r[sp_idx_val]);
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

        if (mapCache.find(devType) == mapCache.end()) {
            auto m = std::make_shared<FieldMap>();
            if (!m->loadFromBinary(ws + "fields_bin/" + binName)) continue;
            mapCache[devType] = m;
        }
        MagneticElement el;
        el.map            = mapCache[devType];
        el.type           = FieldType::MAGNETIC;
        el.useSymmetry    = true;
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
        designPolyline.reserve(els.size() + 2);
        // prelude: первая точка - first_tangent × 10м
        const auto& first = els.front();
        FieldVector t0 = {first.transform.m[0][0],
                          first.transform.m[1][0],
                          first.transform.m[2][0]};
        FieldVector first_pos = {first.transform.tx, first.transform.ty, first.transform.tz};
        designPolyline.push_back({first_pos.x - 10.0*t0.x,
                                  first_pos.y - 10.0*t0.y,
                                  first_pos.z - 10.0*t0.z});
        // тела элементов
        for (const auto& el : els) {
            designPolyline.push_back({el.transform.tx, el.transform.ty, el.transform.tz});
        }
        // postlude: последняя точка + last_tangent × 20м
        const auto& last = els.back();
        FieldVector tN = {last.transform.m[0][0],
                          last.transform.m[1][0],
                          last.transform.m[2][0]};
        FieldVector last_pos = {last.transform.tx, last.transform.ty, last.transform.tz};
        designPolyline.push_back({last_pos.x + 20.0*tN.x,
                                  last_pos.y + 20.0*tN.y,
                                  last_pos.z + 20.0*tN.z});
        std::cout << "  design polyline: " << designPolyline.size() << " pts ("
                  << els.size() << " elements + 2 extensions)\n";
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
    std::cout << "  particles: " << beam.size() << "\n";

    // ---- 6. Magnet positions for loss ----
    std::set<int> spFilter(cfg.lossSpNumbers.begin(), cfg.lossSpNumbers.end());
    std::set<std::string> devExclude(cfg.lossExcludeDevices.begin(),
                                     cfg.lossExcludeDevices.end());
    std::vector<FieldVector> magnetPos;
    std::vector<int>         magnetSp;
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
        double x = std::stod(r[lay_x]);
        double y = std::stod(r[lay_y]);
        double z = std::stod(r[lay_z]);
        auto key = std::make_tuple(std::round(x*1e6)/1e6,
                                   std::round(y*1e6)/1e6,
                                   std::round(z*1e6)/1e6);
        if (seenPos.count(key)) continue;
        seenPos.insert(key);
        magnetPos.push_back({x, y, z});
        magnetSp.push_back(spn);
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
        magnetSp.push_back(-1);  // -1 = explicit target (не из sp_number-фильтра)
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
    auto runTrackingAndLoss = [&](double& outLoss, double& outMeanMM, double& outMaxMM) -> bool {
        for (size_t i = 0; i < beam.size(); ++i) beam[i] = beam0[i];
        std::vector<FinishPlane> planes(beam.size());
        for (size_t i = 0; i < beam.size(); ++i)
            planes[i].setup(fp1, fp2, fp3, startRef, maxTurns);

        const size_t M = magnetPos.size();
        const double INF = std::numeric_limits<double>::infinity();
        std::vector<std::vector<double>> minD2(beam.size(), std::vector<double>(M, INF));
        std::vector<char> reached(beam.size(), 0);

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
                    if (d2 < md[m]) md[m] = d2;
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
        for (const auto& md : minD2) {
            for (size_t m = 0; m < md.size(); ++m) {
                double v = md[m];
                if (!std::isfinite(v)) return false;
                double d = std::sqrt(v);
                double w = magnetWeight[m];
                weighted_sum_d2 += v * w;
                weighted_sum_d  += d * w;
                sum_d += d;
                if (d > max_d) max_d = d;
                total++;
            }
        }
        outLoss   = cfg.lossSquared ? weighted_sum_d2 : weighted_sum_d;
        outMeanMM = (total > 0) ? (sum_d / total) * 1000.0 : 0.0;
        outMaxMM  = max_d * 1000.0;

        // ============== УГЛОВАЯ КОМПОНЕНТА LOSS ==============
        // Σ_по_частицам sin²(θ_i) × dirWeight, где θ_i -- угол между
        // финальным импульсом частицы и эталонным направлением.
        if (cfg.dirWeight > 0.0) {
            double dir_sum = 0.0;
            for (size_t i = 0; i < beam.size(); ++i) {
                const auto& p = beam[i].mom;
                double pn = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
                if (pn < 1e-30) continue;
                double dot = (p.x * cfg.dirReference.x +
                              p.y * cfg.dirReference.y +
                              p.z * cfg.dirReference.z) / pn;
                // Clamp на случай рукости плавающей арифметики
                if (dot > 1.0) dot = 1.0;
                if (dot < -1.0) dot = -1.0;
                double sin2 = 1.0 - dot * dot;
                dir_sum += sin2;
            }
            outLoss += cfg.dirWeight * dir_sum;
        }

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
        if (cfg.envWeight > 0.0 && designPolyline.size() >= 2) {
            double env_sum = 0.0;
            for (size_t i = 0; i < beam.size(); ++i) {
                double a_init  = distanceToPolyline(beam0[i].pos, designPolyline);
                double a_final = distanceToPolyline(beam[i].pos,  designPolyline);
                double growth = a_final - a_init;
                if (growth > 0.0) env_sum += growth * growth;
            }
            outLoss += cfg.envWeight * env_sum;
        }

        return true;
    };

    // Старая runTracking -- используется ТОЛЬКО для финального вывода .bin
    // (в нём всё ещё хранится полная траектория, чтобы visualize.py мог отрисовать).
    auto runTracking = [&](std::vector<std::vector<FieldVector>>& results) -> bool {
        // reset particles
        for (size_t i = 0; i < beam.size(); ++i) beam[i] = beam0[i];
        results.assign(beam.size(), {});
        FinishPlane fPlane;
        // Поскольку FinishPlane хранит callCount per-instance, а мы запускаем
        // 5 частиц параллельно, делаем 5 копий.
        std::vector<FinishPlane> planes(beam.size());
        for (size_t i = 0; i < beam.size(); ++i) {
            planes[i].setup(fp1, fp2, fp3, startRef, maxTurns);
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
    if (!runTrackingAndLoss(curLoss, curMeanMM, curMaxMM)) {
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
        if (!runTrackingAndLoss(newLoss, meanMM, maxMM)) {
            applyGroupValue(rg, oldVal);
            failed++;
            HistEntry e{it, rg.cfg.name, oldVal, newVal, curLoss, bestLoss, curMeanMM, curMaxMM, "FAIL"};
            history.push_back(e);
            std::cout << "[" << it << "/" << cfg.nIter << "] " << rg.cfg.name
                      << "  " << oldVal << " -> " << newVal << "  FAIL  best=" << bestLoss << "\n";
            continue;
        }
        std::string mark;
        if (newLoss < curLoss) {
            curLoss = newLoss;
            curMeanMM = meanMM; curMaxMM = maxMM;
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
    double loopSec = loopTimer.sec();

    // ---- 13. Restore best state ----
    for (size_t k = 0; k < rgroups.size(); ++k) applyGroupValue(rgroups[k], bestValues[k]);

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
