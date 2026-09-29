// matrix_ring.h — АВТОНОМНЫЙ слой кольца матричного ядра (C++): чтение БД, глобальная геометрия,
// инжекция пучка, стыковочные дрейфы, цикл оборотов + ВЧ. Поверх физдвижка matrix_core.h.
//
// Дословный порт геометрии Python-эталона (matrix_tracker.frame/_rot/build_line/check_junctions/
// modified_drift + beam_injection.inject). xref/ref_bow НЕ нужны — их использует только box-тело,
// которого в секторном режиме нет. Валидация — test_matrix_core_ring.py (бит-в-бит с Python).
#pragma once
#include "matrix_core.h"
#include <vector>
#include <string>
#include <map>
#include <array>
#include <fstream>
#include <sstream>
#include <algorithm>

namespace mc {

// ---- 3D вектор (операции как numpy: dot=@, cross, norm) ----------------------------------------
struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() {}
    Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
};
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline Vec3 normalize(const Vec3& a) { double n = norm(a); return (n > 0) ? a * (1.0 / n) : a; }

// Поворот вектора v вокруг ЕДИНИЧНОЙ оси n на угол a (Родриг) — matrix_tracker._rot.
inline Vec3 rot(const Vec3& v, const Vec3& n, double a) {
    double c = std::cos(a), s = std::sin(a);
    return v * c + cross(n, v) * s + n * (dot(n, v) * (1.0 - c));
}
// Ортонормальная тройка (axis,u,v) из углов установки (Tait-Bryan yaw→pitch→roll) — matrix_tracker.frame.
inline void frame(double yaw, double pitch, double roll, Vec3& axis, Vec3& u, Vec3& v) {
    axis = {1, 0, 0}; u = {0, 1, 0}; v = {0, 0, 1};
    if (yaw)   { Vec3 a2 = rot(axis, v, yaw),   u2 = rot(u, v, yaw);   axis = a2; u = u2; }
    if (pitch) { Vec3 a2 = rot(axis, u, pitch), v2 = rot(v, u, pitch); axis = a2; v = v2; }
    if (roll)  { Vec3 u2 = rot(u, axis, roll),  v2 = rot(v, axis, roll); u = u2; v = v2; }
}

struct Frame { Vec3 p, axis, u, v; };

// Параметры ядра прогона (matrix_tracker.Settings): дефолты = matrix_tracker.Settings.DEFAULTS.
struct Settings {
    double loss_radius_m     = 0.012;
    double plane_same_tol_m  = 1e-4;    // плоскости ближе — одна (снап, без пролёта)
    double plane_parallel_tol = 1e-9;   // |axis_A × axis_B| ниже — нормали параллельны
    int    n_slices          = 40;      // фолбэк box-модели (сектор/кик всегда 5)
};

// Узел кольца: физэлемент (тело) + глобальная геометрия. Плоскость инжекции — тоже узел (is_plane, тело=id).
struct Node {
    std::string dev, kind;
    bool   is_plane = false, monitor = false, sector = false, same_plane_as_prev = false;
    Element phys;
    Vec3   center, axis, uhat, vhat, entry, exit;   // «ящик» (для стыков)
    Frame  fr_in, fr_out;                            // рамки стыковки (у сектора — криволинейные)
    // Плоскость инжекции — тело=тождество; иначе физдвижок (Element::body сам даёт дрейф для не-sector/не-kick).
    State  body(const State& s, IsrRng* rng = nullptr) const { return is_plane ? s : phys.body(s, rng); }
};

// ---- CSV (';' разделитель, EU-запятая→точка, срез кавычек/пробелов) — как src/main.cpp ---------
struct Table {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    int col(const std::string& name) const {
        for (size_t i = 0; i < header.size(); ++i) if (header[i] == name) return (int)i;
        return -1;
    }
};
inline std::vector<std::string> splitCSV(const std::string& line) {
    std::vector<std::string> row; std::stringstream ss(line); std::string val;
    while (std::getline(ss, val, ';')) {
        val.erase(0, val.find_first_not_of(" \t\r\n\""));
        size_t e = val.find_last_not_of(" \t\r\n\"");
        val = (e == std::string::npos) ? "" : val.substr(0, e + 1);
        for (char& c : val) if (c == ',') c = '.';
        row.push_back(val);
    }
    return row;
}
inline Table readTable(const std::string& path) {
    Table t; std::ifstream f(path); std::string line;
    if (std::getline(f, line)) t.header = splitCSV(line);
    while (std::getline(f, line)) { if (!line.empty()) t.rows.push_back(splitCSV(line)); }
    return t;
}
inline std::string cell(const std::vector<std::string>& r, int c) {
    return (c >= 0 && c < (int)r.size()) ? r[c] : std::string();
}
inline double cellD(const std::vector<std::string>& r, int c, double def = 0.0) {
    std::string s = cell(r, c);
    if (s.empty() || s == "-") return def;
    try { return std::stod(s); } catch (...) { return def; }
}

// ---- ПЛОСКОСТЬ ИНЖЕКЦИИ (PlaneMagnet) ---------------------------------------------------------
inline Node makePlane(const Vec3& center, const Vec3& axisIn) {
    Node n; n.is_plane = true; n.dev = "INJ_PLANE"; n.kind = "";
    Vec3 a = normalize(axisIn);
    Vec3 u{ -a.y, a.x, 0.0 };
    double nu = norm(u); u = (nu > 1e-9) ? u * (1.0 / nu) : Vec3{0, 1, 0};
    Vec3 v = cross(a, u);
    n.center = center; n.axis = a; n.uhat = u; n.vhat = v; n.entry = center; n.exit = center;
    n.fr_in = { center, a, u, v }; n.fr_out = { center, a, u, v };
    return n;
}

// ---- СБОРКА КОЛЬЦА из БД (build_line + check_junctions) ----------------------------------------
// nsl/order — срезы и порядок сектор/кик-тела (по умолчанию как elegant: 5 и 4); isr — квантовое возбуждение.
inline std::vector<Node> buildLine(const std::string& ws, int lattice, double design_mev, bool sr,
                                   const Settings& stg, int nsl_sector = SECTOR_NSLICES, int yoshida_order = 4,
                                   bool isr = false) {
    Table lc = readTable(ws + "tables/lattice_configs.csv");
    Table cn = readTable(ws + "tables/element_connectivity.csv");
    int cLat = lc.col("lattice_id"), cInst = lc.col("instance_id");
    int cx = lc.col("x"), cy = lc.col("y"), cz = lc.col("z"), cyaw = lc.col("yaw"),
        cpit = lc.col("pitch"), crol = lc.col("roll"), cdev = lc.col("device_type_id"),
        cB = lc.col("mat_B"), cG = lc.col("mat_G"), cGv = lc.col("mat_Gv"), cS = lc.col("mat_S"),
        cO = lc.col("mat_O"), cLeff = lc.col("eff_length"), cKind = lc.col("kind"),
        cNsl = lc.col("n_slices"), cE1 = lc.col("edge1"), cE2 = lc.col("edge2");
    // индекс lattice_configs по instance_id (в пределах решётки)
    std::map<std::string, const std::vector<std::string>*> conf;
    for (auto& r : lc.rows) if (cell(r, cLat) == std::to_string(lattice)) conf[cell(r, cInst)] = &r;
    // порядок из element_connectivity (сортировка по instance_id как int)
    int nLat = cn.col("lattice_id"), nInst = cn.col("instance_id");
    std::vector<std::string> order;
    for (auto& r : cn.rows) if (cell(r, nLat) == std::to_string(lattice)) order.push_back(cell(r, nInst));
    std::sort(order.begin(), order.end(), [](const std::string& a, const std::string& b) {
        return std::stoi(a) < std::stoi(b);
    });
    std::vector<Node> mags;
    size_t i = 0;
    while (i < order.size()) {
        const auto* r0 = conf[order[i]];
        if (!r0) { ++i; continue; }
        double B = cellD(*r0, cB), G = cellD(*r0, cG);
        double Gv = cell(*r0, cGv).empty() ? cellD(*r0, cG) : cellD(*r0, cGv);
        double S = cellD(*r0, cS), O = cellD(*r0, cO);
        std::string kx = cell(*r0, cx), ky = cell(*r0, cy), kz = cell(*r0, cz);   // ключ со-локации (СТРОКИ)
        while (i + 1 < order.size()) {                                            // суммируем со-локованные обмотки
            const auto* rn = conf[order[i + 1]];
            if (rn && cell(*rn, cx) == kx && cell(*rn, cy) == ky && cell(*rn, cz) == kz) {
                B += cellD(*rn, cB); G += cellD(*rn, cG);
                Gv += cell(*rn, cGv).empty() ? cellD(*rn, cG) : cellD(*rn, cGv);
                S += cellD(*rn, cS); O += cellD(*rn, cO);
                ++i;
            } else break;
        }
        Node nd;
        nd.dev = cell(*r0, cdev); nd.kind = cell(*r0, cKind);
        double Leff = cellD(*r0, cLeff);
        int nsl = nsl_sector;       // сектор/кик-тело (_nsl_sector); n_slices из БД — для box-тела (кик корректора)
        double e1 = cellD(*r0, cE1), e2 = cellD(*r0, cE2);
        int nsl_box = (int)cellD(*r0, cNsl, 40.0);
        nd.phys = Element::make(nd.kind, Leff, B, G, Gv, S, O, e1, e2, nsl, design_mev, sr, yoshida_order, isr,
                                nsl_box > 0 ? nsl_box : 40);
        nd.sector = nd.phys.sector; nd.monitor = (nd.kind == "monitor");
        double yaw = cellD(*r0, cyaw), pitch = cellD(*r0, cpit), roll = cellD(*r0, crol);
        nd.center = { cellD(*r0, cx), cellD(*r0, cy), cellD(*r0, cz) };
        frame(yaw, pitch, roll, nd.axis, nd.uhat, nd.vhat);
        nd.entry = nd.center - nd.axis * (0.5 * Leff);
        nd.exit  = nd.center + nd.axis * (0.5 * Leff);
        if (nd.sector) {                                            // криволинейные рамки (⊥ касательной дуги)
            double rho = Leff / nd.phys.theta0;
            Vec3 Oc = nd.center + nd.uhat * rho;                    // центр кривизны
            double th2 = nd.phys.theta0 / 2.0;
            Vec3 entry_pt = Oc + rot(nd.center - Oc, nd.vhat, -th2);
            Vec3 exit_pt  = Oc + rot(nd.center - Oc, nd.vhat,  th2);
            nd.fr_in  = { entry_pt, rot(nd.axis, nd.vhat, -th2), rot(nd.uhat, nd.vhat, -th2), nd.vhat };
            nd.fr_out = { exit_pt,  rot(nd.axis, nd.vhat,  th2), rot(nd.uhat, nd.vhat,  th2), nd.vhat };
        } else {
            nd.fr_in  = { nd.entry, nd.axis, nd.uhat, nd.vhat };
            nd.fr_out = { nd.exit,  nd.axis, nd.uhat, nd.vhat };
        }
        mags.push_back(nd); ++i;
    }
    // check_junctions: снап совпадающих плоскостей (same_plane_as_prev) — по «ящику» (axis/entry/exit).
    for (size_t j = 1; j < mags.size(); ++j) {
        Vec3 nA = mags[j - 1].axis, nB = mags[j].axis, pA = mags[j - 1].exit, pB = mags[j].entry;
        double gap = std::fabs(dot(nB, pA - pB));
        double s = norm(cross(nA, nB));
        mags[j].same_plane_as_prev = (s < stg.plane_parallel_tol) && (gap < stg.plane_same_tol_m);
    }
    return mags;
}

// ---- СТЫКОВОЧНЫЙ ДРЕЙФ (modified_drift) --------------------------------------------------------
constexpr double PARA_TOL = 1e-4;   // |axis_A × axis_B| ниже ⇒ оси параллельны (matrix_tracker._PARA_TOL)
inline State modified_drift(const State& st, const Frame& A, const Frame& B, bool same_plane) {
    Vec3 aP = A.p, aAx = A.axis, aU = A.u, aV = A.v;
    Vec3 bP = B.p, bAx = B.axis, bU = B.u, bV = B.v;
    double d = st.delta;
    Vec3 P = aP + aU * st.x + aV * st.y;
    double ps = std::sqrt((1.0 + d) * (1.0 + d) - st.px * st.px - st.py * st.py);
    Vec3 dir = (aU * st.px + aV * st.py + aAx * ps) * (1.0 / (1.0 + d));
    Vec3 Q;
    if (same_plane) {
        Q = P;
    } else if (norm(cross(aAx, bAx)) < PARA_TOL) {                  // оси параллельны ⇒ чистый дрейф
        Q = P + dir * (dot(bP - P, aAx) / dot(dir, aAx));
    } else {
        double denom = dot(bAx, dir);                              // клин (диполь → элемент под углом)
        if (denom <= 1e-9) throw ParticleLost("beam won't reach next plane");
        double tt = dot(bAx, bP - P) / denom;
        if (tt < -1e-4) throw ParticleLost("next entry plane behind particle");
        Q = P + dir * tt;
    }
    double tt = dot(Q - P, dir);
    Vec3 rel = Q - bP;
    return { dot(rel, bU), (1.0 + d) * dot(dir, bU), dot(rel, bV), (1.0 + d) * dot(dir, bV), st.t + tt, d };
}

// ---- ПУЧОК: загрузка группы, плоскость, место инжекции, канонические координаты ----------------
struct Particle { Vec3 pos, vel; double m, q; };

inline std::vector<Particle> loadGroup(const std::string& ws, int group) {
    Table pg = readTable(ws + "tables/particle_groups.csv");
    int g = pg.col("group_id"), cx = pg.col("x"), cy = pg.col("y"), cz = pg.col("z"),
        vx = pg.col("vx"), vy = pg.col("vy"), vz = pg.col("vz"), cm = pg.col("m"), cq = pg.col("q");
    std::vector<Particle> out;
    for (auto& r : pg.rows) if (cell(r, g) == std::to_string(group)) {
        out.push_back({ { cellD(r, cx), cellD(r, cy), cellD(r, cz) },
                        { cellD(r, vx), cellD(r, vy), cellD(r, vz) }, cellD(r, cm), cellD(r, cq) });
    }
    return out;
}

inline double gamma_of(double v) { return 1.0 / std::sqrt(1.0 - (v / C_LIGHT) * (v / C_LIGHT)); }
// опорный импульс P0 = bm.P_MOM = γ₀·m_e·v₀ (сурогат) — как bip_matrix.
inline double p0_ref() { return gamma_of(V_BEAM) * ME * V_BEAM; }

inline State toCanonical(const Particle& p, const Node& plane) {
    Vec3 rel = p.pos - plane.center;
    double x = dot(rel, plane.uhat), y = dot(rel, plane.vhat);
    double v = norm(p.vel); Vec3 dir = p.vel * (1.0 / v);
    double delta = (gamma_of(v) * p.m * v - p0_ref()) / p0_ref();
    double px = (1.0 + delta) * dot(dir, plane.uhat), py = (1.0 + delta) * dot(dir, plane.vhat);
    return { x, px, y, py, 0.0, delta };
}

// Результат инжекции: плоскость-старт, упорядоченное кольцо, финиш, канон. состояния, замкнутость, периметр.
struct Injection { Node plane; std::vector<Node> ordered; Node finish; std::vector<State> states;
                   bool closed; double ring_len = 0.0; std::vector<double> entry_arc; };

inline Injection inject(const std::string& ws, int lattice, int group, double design_mev, bool sr,
                        const Settings& stg) {
    std::vector<Node> mags = buildLine(ws, lattice, design_mev, sr, stg);
    std::vector<Particle> parts = loadGroup(ws, group);
    if (parts.empty()) throw std::runtime_error("beam group not found");
    // плоскость пучка: центроид позиций + средняя нормированная скорость
    Vec3 csum, vsum, vdir;
    for (auto& p : parts) { csum = csum + p.pos; vsum = vsum + p.vel; vdir = vdir + normalize(p.vel); }
    Vec3 centroid = csum * (1.0 / parts.size());
    Vec3 normal = normalize(vsum * (1.0 / parts.size()));
    vdir = normalize(vdir * (1.0 / parts.size()));
    // поиск дрейфа, где орбита пересекает плоскость пучка по направлению скорости
    int N = (int)mags.size(); int bi = -1; Vec3 bcenter, baxis; double bdist = 0;
    for (int k = 0; k < N; ++k) {
        Vec3 A = mags[k].exit, B = mags[(k + 1) % N].entry, AB = B - A;
        double L = norm(AB);
        if (L < 1e-9) continue;
        Vec3 axis = AB * (1.0 / L);
        double denom = dot(normal, AB);
        if (std::fabs(denom) < 1e-12) continue;
        double t = dot(normal, centroid - A) / denom;
        if (t >= -0.02 && t <= 1.02) {
            Vec3 center = A + AB * t;
            if (dot(vdir, axis) <= 0) continue;
            double dist = norm(center - centroid);
            if (bi < 0 || dist < bdist) { bi = k; bcenter = center; baxis = axis; bdist = dist; }
        }
    }
    Injection inj;
    if (bi >= 0) {                                             // ЗАМКНУТОЕ кольцо
        inj.plane = makePlane(bcenter, baxis);
        for (int k = bi + 1; k < N; ++k) inj.ordered.push_back(mags[k]);
        for (int k = 0; k <= bi; ++k) inj.ordered.push_back(mags[k]);
        inj.finish = inj.plane; inj.closed = true;
    } else {                                                   // ОТКРЫТЫЙ сегмент
        inj.plane = makePlane(centroid, vdir);
        inj.ordered = mags;
        inj.finish = makePlane(mags.back().exit, mags.back().axis);
        inj.closed = false;
    }
    for (auto& p : parts) inj.states.push_back(toCanonical(p, inj.plane));
    // опорная арк-длина входа каждого элемента (entry_arc) + периметр (ring_len) — как local_trajectories.
    const auto& O = inj.ordered;
    inj.entry_arc.resize(O.size());
    inj.entry_arc[0] = dot(O[0].fr_in.p - inj.plane.fr_out.p, inj.plane.fr_out.axis);   // gap_in
    for (size_t j = 1; j < O.size(); ++j) {
        double ord_gap = dot(O[j].fr_in.p - O[j - 1].fr_out.p, O[j - 1].fr_out.axis);
        inj.entry_arc[j] = inj.entry_arc[j - 1] + O[j - 1].phys.Leff + ord_gap;
    }
    double gap_out = inj.closed ? dot(inj.finish.fr_in.p - O.back().fr_out.p, O.back().fr_out.axis) : 0.0;
    inj.ring_len = inj.entry_arc.back() + O.back().phys.Leff + gap_out;
    return inj;
}

// ---- ТРЕК: N оборотов (плоскость → кольцо → финиш) + ВЧ (ideal: восстановить δ за оборот) -------
struct TrackResult { State final_state; int surv_turns; bool lost; };

inline TrackResult trackParticle(const Injection& inj, const State& st0, int n_turns, bool rf_ideal) {
    State s = st0; int done = 0;
    const auto& O = inj.ordered;
    try {
        for (int turn = 0; turn < n_turns; ++turn) {
            double d0 = s.delta;
            s = modified_drift(s, inj.plane.fr_out, O[0].fr_in, O[0].same_plane_as_prev);
            s = O[0].body(s);
            for (size_t j = 1; j < O.size(); ++j) {
                s = modified_drift(s, O[j - 1].fr_out, O[j].fr_in, O[j].same_plane_as_prev);
                s = O[j].body(s);
            }
            s = modified_drift(s, O.back().fr_out, inj.finish.fr_in, inj.finish.same_plane_as_prev);
            if (rf_ideal) s.delta = d0;
            done = turn + 1;
        }
    } catch (const ParticleLost&) {
        return { s, done, true };
    }
    return { s, done, false };
}

// ДЛИННЫЙ трек: состояние (x,px,y,py) у плоскости старта РАЗ ЗА ОБОРОТ (как long_trajectories: tt[0]=старт,
// tt[k]=после k-го оборота, ≤n_turns+1 точек). Для вывода в .npz анализатору.
struct LongResult { std::vector<std::array<double, 4>> tt; int surv_turns = 0; State final_state; bool lost = false; };

inline LongResult trackLong(const Injection& inj, const State& st0, int n_turns, bool rf_ideal) {
    LongResult r; State s = st0;
    r.tt.push_back({ s.x, s.px, s.y, s.py });
    const auto& O = inj.ordered;
    try {
        for (int turn = 0; turn < n_turns; ++turn) {
            double d0 = s.delta;
            s = modified_drift(s, inj.plane.fr_out, O[0].fr_in, O[0].same_plane_as_prev);
            s = O[0].body(s);
            for (size_t j = 1; j < O.size(); ++j) {
                s = modified_drift(s, O[j - 1].fr_out, O[j].fr_in, O[j].same_plane_as_prev);
                s = O[j].body(s);
            }
            s = modified_drift(s, O.back().fr_out, inj.finish.fr_in, inj.finish.same_plane_as_prev);
            if (rf_ideal) s.delta = d0;
            r.tt.push_back({ s.x, s.px, s.y, s.py });
            r.surv_turns = turn + 1;
        }
    } catch (const ParticleLost&) { r.lost = true; }
    r.final_state = s;
    return r;
}

// Локальные (x,y) в рамке fr → глобальные XYZ (matrix_traj_export.gp).
inline void gp(const Frame& fr, double x, double y, double& X, double& Y, double& Z) {
    Vec3 p = fr.p + fr.u * x + fr.v * y; X = p.x; Y = p.y; Z = p.z;
}

// ПОЛНЫЙ трек: запись (S,x,px,y,py,ct,delta,elem)+глоб.XYZ на ГРАНИЦАХ каждого элемента за все обороты
// (для анализа x(S)/спектра/тюнов) с равномерным прореживанием до target точек. БЕЗ ВЧ — как local_trajectories.
struct FullResult {
    std::vector<double> S, x, px, y, py, ct, delta, X, Y, Z;
    std::vector<int64_t> elem;
    int surv_turns = 0;
    size_t size() const { return S.size(); }
};

inline FullResult trackFull(const Injection& inj, const State& st0, int n_turns, int target) {
    FullResult r; State st = st0;
    const auto& O = inj.ordered; const auto& ea = inj.entry_arc;
    auto rec = [&](double s, const State& p, int64_t el, const Frame& fr) {
        r.S.push_back(s); r.x.push_back(p.x); r.px.push_back(p.px); r.y.push_back(p.y);
        r.py.push_back(p.py); r.ct.push_back(p.t); r.delta.push_back(p.delta); r.elem.push_back(el);
        double X, Y, Z; gp(fr, p.x, p.y, X, Y, Z); r.X.push_back(X); r.Y.push_back(Y); r.Z.push_back(Z);
    };
    rec(0.0, st, 0, inj.plane.fr_out);                                   // старт у плоскости
    try {
        for (int turn = 0; turn < n_turns; ++turn) {
            double base = turn * inj.ring_len;
            st = modified_drift(st, inj.plane.fr_out, O[0].fr_in, O[0].same_plane_as_prev);
            rec(base + ea[0], st, 1, O[0].fr_in);                        // вход эл.0
            st = O[0].body(st);
            rec(base + ea[0] + O[0].phys.Leff, st, 1, O[0].fr_out);      // выход эл.0
            for (size_t j = 1; j < O.size(); ++j) {
                st = modified_drift(st, O[j - 1].fr_out, O[j].fr_in, O[j].same_plane_as_prev);
                rec(base + ea[j], st, (int64_t)(j + 1), O[j].fr_in);
                st = O[j].body(st);
                rec(base + ea[j] + O[j].phys.Leff, st, (int64_t)(j + 1), O[j].fr_out);
            }
            if (inj.closed) {                                           // финиш = плоскость (замыкание)
                st = modified_drift(st, O.back().fr_out, inj.finish.fr_in, inj.finish.same_plane_as_prev);
                rec(base + inj.ring_len, st, 0, inj.finish.fr_in);
            } else {                                                    // ОТКРЫТЫЙ (демо): симметричный ЗАВЕРШАЮЩИЙ дрейф = ведущему ea[0]
                double trail = ea[0];                                   // inject ставит финиш на выходе элемента; отодвигаем на ea[0] (=gap_out при ENTRY=EXIT)
                Frame ftr = inj.finish.fr_in;                           // финиш-плоскость по BOX-оси (как Python bi.inject), НЕ tilted fr_out
                ftr.p = ftr.p + ftr.axis * trail;
                st = modified_drift(st, O.back().fr_out, ftr, false);
                rec(base + inj.ring_len + trail, st, (int64_t)(O.size() + 1), ftr);
            }
            r.surv_turns = turn + 1;
        }
    } catch (const ParticleLost&) {}
    size_t n = r.size();                                                 // равномерное прореживание (шаг ceil(n/target))
    if (target > 0 && (int)n > target) {
        size_t step = (n + (size_t)target - 1) / (size_t)target;
        FullResult ds; ds.surv_turns = r.surv_turns;
        for (size_t i = 0; i < n; i += step) {
            ds.S.push_back(r.S[i]); ds.x.push_back(r.x[i]); ds.px.push_back(r.px[i]); ds.y.push_back(r.y[i]);
            ds.py.push_back(r.py[i]); ds.ct.push_back(r.ct[i]); ds.delta.push_back(r.delta[i]);
            ds.elem.push_back(r.elem[i]); ds.X.push_back(r.X[i]); ds.Y.push_back(r.Y[i]); ds.Z.push_back(r.Z[i]);
        }
        return ds;
    }
    return r;
}

} // namespace mc
