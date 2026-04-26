#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <omp.h>
#include <map>
#include <fstream>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
#endif

#include "AnsysConverter.h"
#include "FieldMap.h"
#include "LatticeManager.h"
#include "Integrator.h"

struct StepTimer {
    std::chrono::high_resolution_clock::time_point start;
    StepTimer() : start(std::chrono::high_resolution_clock::now()) {}
    long long ms() {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    }
};

struct FinishPlane {
    double x_limit;
    void setup(double x) { x_limit = x; }
    bool hasCrossed(const FieldVector& pos) const { return (pos.x >= x_limit); }
};

std::vector<std::vector<std::string>> readCSV(std::string path) {
    std::vector<std::vector<std::string>> data;
    std::ifstream f(path);
    if (!f.is_open()) return data;
    std::string line;
    std::getline(f, line);
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row;
        for (char &c : line) if (c == ';') c = ',';
        std::stringstream ss(line);
        std::string val;
        while (std::getline(ss, val, ',')) {
            val.erase(0, val.find_first_not_of(" \t\r\n\""));
            val.erase(val.find_last_not_of(" \t\r\n\"") + 1);
            row.push_back(val);
        }
        if (!row.empty()) data.push_back(row);
    }
    return data;
}

int main(int argc, char* argv[]) {
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    #endif

    if (argc < 2) {
        std::cout << "Usage: SynchrotronTracker.exe <exp_id>" << std::endl;
        return 1;
    }

    int targetExpId = std::stoi(argv[1]);
    std::string ws = "../workspace/"; 
    auto total_timer = StepTimer();

    // 1. Experiment loading
    auto expData = readCSV(ws + "tables/experiments.csv");
    int latId = -1, grpId = -1, maxTurns = 0;
    double dt = 0, finishX = 5.0;

    for (auto& r : expData) {
        if (std::stoi(r[0]) == targetExpId) {
            latId = std::stoi(r[1]);
            grpId = std::stoi(r[2]);
            dt = std::stod(r[3]);
            maxTurns = std::stoi(r[4]);
            if (r.size() > 5) finishX = std::stod(r[5]);
            break;
        }
    }
    
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
            if (mapCache.find(type) == mapCache.end()) {
                for (auto& h : hwData) {
                    if (h[0] == type) {
                        auto m = std::make_shared<FieldMap>();
                        if (m->loadFromBinary(ws + "fields_bin/" + h[3])) {
                            mapCache[type] = m;
                        }
                        break;
                    }
                }
            }
            if (mapCache.count(type)) {
                MagneticElement el;
                el.map = mapCache[type];
                el.setOrientation(std::stod(r[3]), std::stod(r[4]), std::stod(r[5]), 
                                  std::stod(r[6]), std::stod(r[7]), std::stod(r[8]));
                el.currentScale = std::stod(r[9]) / el.map->getNominalArg();
                el.useSymmetry = true;
                lattice.addElement(el);
            }
        }
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
    FinishPlane fPlane; fPlane.setup(finishX);

    for (int turn = 1; turn <= maxTurns; ++turn) {
        auto turn_timer = StepTimer();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < (int)beam.size(); ++i) {
            int steps = 0;
            while (steps < 2000000) { 
                Integrator::rk4Step(beam[i], masses[i], charges[i], lattice, dt);
                if (steps % 10 == 0) results[i].push_back(beam[i].pos); 
                steps++;
                if (fPlane.hasCrossed(beam[i].pos)) break;
            }
        }
        std::cout << "Turn " << turn << ": " << turn_timer.ms() << " ms" << std::endl;
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
