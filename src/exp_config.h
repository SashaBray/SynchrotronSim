// exp_config.h — конфиг опыта для диффурного ядра (Boris: SynchrotronTracker, SynchrotronOptimizer).
//
// Читает experiments.csv и particle_groups.csv ПО ИМЕНАМ столбцов (раньше — по номерам). Порт
// beam_injection (Python): старт-плоскость start_x..start_nz (точка + нормаль), финиш finish_x..finish_nz
// (пусто ⇒ = старт, кольцо), пучок — локальные (x, xp, y, yp, ct, delta) в рамке старт-плоскости →
// глобальные pos/vel: pos = c + x·u + y·v, направление ∝ xp·u + yp·v + n, |p| = (1+δ)·P0 (P0 = γ₀·m_e·v₀,
// опорный импульс сурогатного пучка — как matrix_core.h/bip_matrix). ct диффурное ядро не использует (ВЧ нет).
// Рамка плоскости — как PlaneMagnet / mc::makePlane: u — горизонталь (n повёрнута на +90° в XY), v = n×u.
#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include "FieldMap.h"   // FieldVector

namespace expcfg {

constexpr double C_LIGHT = 299792458.0;
constexpr double ME      = 9.1093837e-31;
constexpr double V_BEAM  = 2.967e8;

struct Table {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    int col(const std::string& name) const {
        for (size_t i = 0; i < header.size(); ++i) if (header[i] == name) return (int)i;
        return -1;
    }
    std::string cell(const std::vector<std::string>& r, const std::string& name) const {
        int c = col(name);
        return (c >= 0 && c < (int)r.size()) ? r[c] : std::string();
    }
};

inline std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> row; std::stringstream ss(line); std::string val;
    while (std::getline(ss, val, ';')) {
        val.erase(0, val.find_first_not_of(" \t\r\n\""));
        size_t e = val.find_last_not_of(" \t\r\n\"");
        val = (e == std::string::npos) ? "" : val.substr(0, e + 1);
        for (char& c : val) if (c == ',') c = '.';               // EU-запятая → точка
        row.push_back(val);
    }
    return row;
}

inline Table read(const std::string& path) {
    Table t; std::ifstream f(path); std::string line;
    if (!f.is_open()) throw std::runtime_error("нет файла " + path);
    if (std::getline(f, line)) t.header = split(line);
    while (std::getline(f, line)) if (!line.empty()) t.rows.push_back(split(line));
    return t;
}

inline FieldVector vadd(const FieldVector& a, const FieldVector& b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
inline FieldVector vmul(const FieldVector& a, double s) { return { a.x * s, a.y * s, a.z * s }; }
inline double vdot(const FieldVector& a, const FieldVector& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline FieldVector vcross(const FieldVector& a, const FieldVector& b) {
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}
inline FieldVector vunit(const FieldVector& a) { double n = std::sqrt(vdot(a, a)); return vmul(a, 1.0 / n); }

struct Plane { FieldVector c, n, u, v; };

inline Plane makePlane(const FieldVector& c, const FieldVector& nIn) {
    Plane p; p.c = c; p.n = vunit(nIn);
    FieldVector u = { -p.n.y, p.n.x, 0.0 };
    double nu = std::sqrt(vdot(u, u));
    p.u = (nu > 1e-9) ? vmul(u, 1.0 / nu) : FieldVector{ 0.0, 1.0, 0.0 };
    p.v = vcross(p.n, p.u);
    return p;
}

// Плоскость из 6 столбцов prefix{x,y,z,nx,ny,nz}: все пусты ⇒ false; частично ⇒ ошибка.
inline bool readPlane(const Table& t, const std::vector<std::string>& r, const std::string& pre, Plane& out) {
    static const char* K[6] = { "x", "y", "z", "nx", "ny", "nz" };
    double v[6]; int filled = 0;
    for (int k = 0; k < 6; ++k) {
        std::string s = t.cell(r, pre + K[k]);
        if (!s.empty() && s != "-") { v[k] = std::stod(s); ++filled; } else v[k] = 0.0;
    }
    if (filled == 0) return false;
    if (filled != 6) throw std::runtime_error("плоскость " + pre + ": заполнены не все 6 столбцов");
    out = makePlane({ v[0], v[1], v[2] }, { v[3], v[4], v[5] });
    return true;
}

struct Experiment {
    int lattice = -1, group = -1, maxTurns = 0, trajStride = 1;
    double dt = 0.0, maxTime = 0.0;
    Plane start, finish; bool hasFinish = false;
};

inline Experiment readExperiment(const std::string& ws, int exp) {
    Table t = read(ws + "tables/experiments.csv");
    if (t.col("start_x") < 0)
        throw std::runtime_error("experiments.csv в старом формате (p1..p3/dir) — запустите migrate_explicit_planes.py");
    for (auto& r : t.rows) {
        if (t.cell(r, "exp_id") != std::to_string(exp)) continue;
        Experiment e;
        e.lattice = std::stoi(t.cell(r, "lattice_id"));
        e.group = std::stoi(t.cell(r, "group_id"));
        std::string s;
        if (!(s = t.cell(r, "dt")).empty()) e.dt = std::stod(s);
        if (!(s = t.cell(r, "max_turns")).empty()) e.maxTurns = std::stoi(s);
        if (!(s = t.cell(r, "traj_stride")).empty()) { try { e.trajStride = std::max(1, std::stoi(s)); } catch (...) {} }
        if (!(s = t.cell(r, "max_time_s")).empty() && s != "-") { try { e.maxTime = std::stod(s); } catch (...) {} }
        if (!readPlane(t, r, "start_", e.start))
            throw std::runtime_error("опыт " + std::to_string(exp) + ": не задана старт-плоскость (exp_planes.py)");
        e.hasFinish = readPlane(t, r, "finish_", e.finish);
        if (!e.hasFinish) e.finish = e.start;                      // пусто ⇒ кольцо: финиш = старт
        return e;
    }
    throw std::runtime_error("Experiment " + std::to_string(exp) + " not found");
}

struct Particle { FieldVector pos, vel; double m = 0, q = 0; };

// Пучок группы в рамке старт-плоскости → глобальные pos/vel (beam_injection.local_to_global).
inline std::vector<Particle> readBeam(const std::string& ws, int group, const Plane& start) {
    Table t = read(ws + "tables/particle_groups.csv");
    if (t.col("xp") < 0)
        throw std::runtime_error("particle_groups.csv в старом формате — запустите migrate_explicit_planes.py");
    const double g0 = 1.0 / std::sqrt(1.0 - (V_BEAM / C_LIGHT) * (V_BEAM / C_LIGHT));
    const double P0 = g0 * ME * V_BEAM;
    auto num = [&](const std::vector<std::string>& r, const char* k) {
        std::string s = t.cell(r, k); return (s.empty() || s == "-") ? 0.0 : std::stod(s);
    };
    std::vector<Particle> out;
    for (auto& r : t.rows) {
        if (t.cell(r, "group_id") != std::to_string(group)) continue;
        Particle p; p.m = num(r, "m"); p.q = num(r, "q");
        double x = num(r, "x"), xp = num(r, "xp"), y = num(r, "y"), yp = num(r, "yp"), d = num(r, "delta");
        p.pos = vadd(start.c, vadd(vmul(start.u, x), vmul(start.v, y)));
        FieldVector dir = vunit(vadd(vadd(vmul(start.u, xp), vmul(start.v, yp)), start.n));
        double gb = (1.0 + d) * P0 / (p.m * C_LIGHT);              // γβ
        p.vel = vmul(dir, C_LIGHT * gb / std::sqrt(1.0 + gb * gb));
        out.push_back(p);
    }
    return out;
}

// Финиш-плоскость в виде FinishPlane: три точки c±u+v, c−v и направление пересечения = нормаль.
inline void finishPoints(const Plane& f, FieldVector& p1, FieldVector& p2, FieldVector& p3, FieldVector& dir) {
    p1 = vadd(f.c, vadd(f.u, f.v));
    p2 = vadd(f.c, vadd(vmul(f.u, -1.0), f.v));
    p3 = vadd(f.c, vmul(f.v, -1.0));
    dir = f.n;
}

} // namespace expcfg
