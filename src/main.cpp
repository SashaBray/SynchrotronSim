#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <omp.h>
#include <map>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <system_error>

#ifdef _WIN32
#include <windows.h>
#endif

#include "AnsysConverter.h"
#include "FieldMap.h"
#include "LatticeManager.h"
#include "Integrator.h"
#include "FinishPlane.h"

struct StepTimer {
    std::chrono::high_resolution_clock::time_point start;
    StepTimer() : start(std::chrono::high_resolution_clock::now()) {}
    long long ms() {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    }
};

std::vector<std::vector<std::string>> readCSV(std::string path) {
    std::vector<std::vector<std::string>> data;
    std::ifstream f(path);
    if (!f.is_open()) return data;
    std::string line;
    std::getline(f, line); // Skip header
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string val;
        // CSV использует ';' как разделитель полей, чтобы не путаться с EU-форматом '1,00E-11'.
        while (std::getline(ss, val, ';')) {
            val.erase(0, val.find_first_not_of(" \t\r\n\""));
            val.erase(val.find_last_not_of(" \t\r\n\"") + 1);
            // EU-формат: запятая как десятичный разделитель -> приводим к '.'
            for (char &c : val) if (c == ',') c = '.';
            row.push_back(val);
        }
        if (!row.empty()) data.push_back(row);
    }
    return data;
}

// Поиск единственной папки workspace. Поднимаемся вверх от каталога с .exe
// (а также от argv[0] и текущей рабочей директории), пока не найдём папку,
// содержащую подкаталог "workspace". Это позволяет запускать .exe из любого места.
std::string findWorkspace(const char* argv0) {
    namespace fs = std::filesystem;
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
            if (fs::is_directory(candidate, ec)) {
                return candidate.generic_string() + "/";
            }
            fs::path parent = dir.parent_path();
            if (parent.empty() || parent == dir) break;
            dir = parent;
        }
    }
    return "workspace/"; // запасной вариант
}

int main(int argc, char* argv[]) {
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    #endif

    if (argc < 2) {
        std::cout << "Usage:" << std::endl;
        std::cout << "  SynchrotronTracker.exe <exp_id>   - run tracking" << std::endl;
        std::cout << "  SynchrotronTracker.exe --convert  - regenerate all .bin field maps from .fld" << std::endl;
        return 1;
    }

    std::string ws = findWorkspace(argv[0]);
    std::cout << "Workspace: " << ws << std::endl;

    // Режим конверсии: .fld -> .bin по hardware_library.csv
    if (std::string(argv[1]) == "--convert") {
        auto hwData = readCSV(ws + "tables/hardware_library.csv");
        int converted = 0, failed = 0;
        for (auto& h : hwData) {
            if (h.size() < 5 || h[0].empty()) continue;
            std::string rawPath = ws + "fields_raw/" + h[3];
            std::string binPath = ws + "fields_bin/" + h[4];
            double nominalArg = 1.0;
            if (h.size() > 7 && h[7] != "-") {
                try { nominalArg = std::stod(h[7]); } catch (...) { nominalArg = 1.0; }
            }
            std::cout << "[" << h[0] << "] " << h[3] << " -> " << h[4]
                      << " (nominal=" << nominalArg << ")" << std::endl;
            if (AnsysConverter::convert(rawPath, binPath, nominalArg)) converted++;
            else { failed++; std::cerr << "  FAILED" << std::endl; }
        }
        std::cout << "Done. Converted: " << converted << ", Failed: " << failed << std::endl;
        return 0;
    }

    int targetExpId = std::stoi(argv[1]);
    auto total_timer = StepTimer();

    // 1. Experiment loading
    auto expData = readCSV(ws + "tables/experiments.csv");
    int latId = -1, grpId = -1, maxTurns = 0;
    int trajStride = 1;
    double dt = 0;
    // Финишная плоскость задаётся тремя точками p1, p2, p3
    FieldVector fp1 = {5.0, -1.0, -1.0};
    FieldVector fp2 = {5.0,  1.0, -1.0};
    FieldVector fp3 = {5.0,  0.0,  1.0};

    for (auto& r : expData) {
        if (std::stoi(r[0]) == targetExpId) {
            latId = std::stoi(r[1]);
            grpId = std::stoi(r[2]);
            dt = std::stod(r[3]);
            maxTurns = std::stoi(r[4]);
            // Колонка 5 -- traj_stride: сохранять каждую N-ю точку траектории
            // (проверка финишной плоскости всё равно идёт каждый шаг).
            // 1 = старое поведение, N>1 = прорежение. Backward compat:
            // если в строке всего 14 столбцов -- считаем traj_stride=1, а
            // p1..p3 берём по старому смещению (с 5-го).
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
    std::cout << "  traj_stride=" << trajStride
              << " (save 1 of every " << trajStride << " trajectory points)" << std::endl;
    
    if (latId == -1) {
        std::cerr << "Error: Experiment " << targetExpId << " not found." << std::endl;
        return 1;
    }

    // 2. Lattice setup
    StepTimer init_timer;
    AcceleratorConfig lattice;
    auto latData = readCSV(ws + "tables/lattice_configs.csv");
    auto hwData = readCSV(ws + "tables/hardware_library.csv");
    std::map<std::string, std::shared_ptr<FieldMap>> mapCache;

    for (auto& r : latData) {
        if (std::stoi(r[0]) == latId) {
            std::string type = r[2];
            
            // Находим параметры типа в библиотеке hardware_library
            bool isActive = true;
            double fieldDir = 1.0;
            double corrCoeff = 1.0;
            std::string binName = "";
            bool foundType = false;

            for (auto& h : hwData) {
                if (h[0] == type) {
                    isActive = (h[1] == "True" || h[1] == "true" || h[1] == "1");
                    binName = h[4];
                    // Колонки CSV: 0:id 1:Active 2:field_type 3:raw 4:bin 5:stl
                    //              6:use_symmetry 7:nominal_arg 8:field_direction 9:correct_coeff
                    fieldDir = std::stod(h[8]);
                    corrCoeff = std::stod(h[9]);
                    foundType = true;
                    break;
                }
            }

            if (!foundType) continue;

            if (mapCache.find(type) == mapCache.end()) {
                auto m = std::make_shared<FieldMap>();
                if (m->loadFromBinary(ws + "fields_bin/" + binName)) {
                    mapCache[type] = m;
                }
            }

            if (mapCache.count(type)) {
                MagneticElement el;
                el.map = mapCache[type];
                el.type = FieldType::MAGNETIC; // По умолчанию B
                el.useSymmetry = true;

                // Переносим калибровочные параметры
                el.isActive = isActive;
                el.fieldDirection = fieldDir;
                // Эффективный correct_coeff = hardware × per-instance (колонка 10
                // в lattice_configs.csv, пустая или "1.0" если не задана).
                double perInstanceCoeff = 1.0;
                if (r.size() > 10 && !r[10].empty() && r[10] != "-") {
                    try { perInstanceCoeff = std::stod(r[10]); } catch (...) {}
                }
                el.correctCoeff = corrCoeff * perInstanceCoeff;

                el.setOrientation(std::stod(r[3]), std::stod(r[4]), std::stod(r[5]),
                                  std::stod(r[6]), std::stod(r[7]), std::stod(r[8]));

                // Если элемент активный, парсим ток. Если постоянный магнит ("-"), ставим 1.0
                if (isActive && r[9] != "-") {
                    el.currentScale = std::stod(r[9]) / el.map->getNominalArg();
                } else {
                    el.currentScale = 1.0;
                }

                lattice.addElement(el);
            }
        }
    }

    // Построить spatial index для быстрого поиска магнитов в окрестности
    // текущей позиции частицы. Без этого getTotalFields() гоняет линейный
    // цикл по всем ~2762 элементам каждый шаг интегрирования.
    {
        StepTimer t_idx;
        int totalRegs = lattice.buildSpatialIndex();
        int nbx, nby, nNonEmpty, maxLen; double bsize, avgLen;
        lattice.getGridStats(nbx, nby, bsize, nNonEmpty, maxLen, avgLen);
        std::cout << "Spatial index: " << t_idx.ms() << " ms, "
                  << nbx << "x" << nby << " buckets (size=" << bsize << " m), "
                  << nNonEmpty << " non-empty, "
                  << "avg=" << avgLen << "/bucket, max=" << maxLen
                  << " (" << totalRegs << " regs total)" << std::endl;
    }

    // 3. Beam loading
    std::vector<ParticleState> beam;
    std::vector<double> masses, charges;
    auto partData = readCSV(ws + "tables/particle_groups.csv");
    for (auto& r : partData) {
        if (std::stoi(r[0]) == grpId) {
            masses.push_back(std::stod(r[2]));
            charges.push_back(std::stod(r[3]));
            ParticleState s;
            s.pos = {std::stod(r[4]), std::stod(r[5]), std::stod(r[6])};
            FieldVector v0 = {std::stod(r[7]), std::stod(r[8]), std::stod(r[9])};
            s.mom = Integrator::velocityToMomentum(v0, masses.back());
            beam.push_back(s);
        }
    }
    std::cout << "Init: " << init_timer.ms() << " ms. Particles: " << beam.size() << std::endl;

    // 4. Tracking
    std::vector<std::vector<FieldVector>> results(beam.size());
    FieldVector startRef = beam.empty() ? FieldVector{0.0, 0.0, 0.0} : beam[0].pos;
    // ВНИМАНИЕ: max_turns теперь репокупируется как "skipSteps" --
    // число вызовов hasCrossed, которые игнорируются перед первой проверкой
    // финишной плоскости. Это нужно для плоскостей, поставленных в начале
    // трассы (например на стартовой точке) -- частица пролетает кольцо и
    // только потом плоскость начинает проверяться.
    // Для коротких экспериментов (плоскость впереди) ставьте max_turns = 0 или 1.
    // У каждой частицы СВОЯ FinishPlane (race-free, плюс per-particle skipSteps).
    std::vector<FinishPlane> planes(beam.size());
    for (size_t i = 0; i < beam.size(); ++i) {
        planes[i].setup(fp1, fp2, fp3, startRef, maxTurns);
    }
    if (!planes.empty()) {
        std::cout << "Finish plane: p0=(" << planes[0].p0.x << "," << planes[0].p0.y << "," << planes[0].p0.z
                  << ")  n=(" << planes[0].n.x << "," << planes[0].n.y << "," << planes[0].n.z << ")"
                  << "  skipSteps=" << planes[0].skipSteps
                  << std::endl;
    }

    {
        auto turn_timer = StepTimer();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < (int)beam.size(); ++i) {
            int steps = 0;
            // ВСЕГДА сохраняем стартовую точку, чтобы по bin'у можно было
            // восстановить начало траектории.
            results[i].push_back(beam[i].pos);
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
        std::cout << "Tracking: " << turn_timer.ms() << " ms" << std::endl;
    }

    // 5. Export
    std::string binP = ws + "results/trajectories_exp" + std::to_string(targetExpId) + ".bin";
    std::ofstream out(binP, std::ios::binary);
    for (size_t i = 0; i < results.size(); ++i) {
        uint32_t sz = (uint32_t)results[i].size();
        out.write((char*)&sz, 4);
        out.write((char*)results[i].data(), sz * sizeof(FieldVector));
        std::cout << "Part " << i+1 << ": " << sz << " points" << std::endl;
    }
    out.close();

    std::string csvP = ws + "results/debug_exp" + std::to_string(targetExpId) + ".csv";
    std::ofstream csv(csvP);
    csv << "part_id,x,y,z\n";
    for(size_t i=0; i<results.size(); ++i)
        for(auto& p : results[i]) csv << i+1 << "," << p.x << "," << p.y << "," << p.z << "\n";
    csv.close();

    std::cout << "Total time: " << total_timer.ms() << " ms" << std::endl;

    return 0;
}
