// matrix_ring.h — АВТОНОМНЫЙ слой кольца матричного ядра (C++): чтение БД, глобальная геометрия,
// инжекция пучка, стыковочные дрейфы, цикл оборотов + ВЧ. Поверх физдвижка matrix_core.h.
//
// Дословный порт геометрии Python-эталона (matrix_tracker.frame/_rot/build_line/check_junctions/
// modified_drift + beam_injection: явные старт/финиш-плоскости опыта, пучок в рамке старта, резонаторы). xref/ref_bow НЕ нужны — их использует только box-тело,
// которого в секторном режиме нет. Валидация — test_matrix_core_ring.py (бит-в-бит с Python).
#pragma once
#include "matrix_core.h"
#include <vector>
#include <string>
#include <map>
#include <cstdlib>
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
    double theta0_0 = 0, edge1_0 = 0, edge2_0 = 0;   // исходные угол/грани изгиба (база BSCALE — живой зазор ИД)
    bool   cavity = false;                           // узел-резонатор пути (плоскость, кик делает трекер ВЧ)
    double cav_share = 0, s_c = 0;                   // доля напряжения, опорная дуга от старт-плоскости
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

// Срезы сектор/кик-тела ПО ТИПУ УСТРОЙСТВА: env SS_NSL_DEV="dev=n,dev=n" перекрывает nsl_sector для
// перечисленных device_type_id (например, короткие слабые полюса ондуляторов — 1 срез). Нет переменной —
// всё как было. Python-эталон читает ту же переменную (matrix_tracker._nsl_by_dev) — паритет бит-в-бит.
inline std::map<std::string, int> nslByDev() {
    std::map<std::string, int> m;
    const char* env = std::getenv("SS_NSL_DEV");
    if (!env) return m;
    std::string s(env);
    size_t i = 0;
    while (i <= s.size()) {
        size_t j = s.find(',', i); if (j == std::string::npos) j = s.size();
        std::string tok = s.substr(i, j - i);
        size_t e = tok.find('=');
        if (e != std::string::npos) {
            int n = std::atoi(tok.c_str() + e + 1);
            if (n > 0) m[tok.substr(0, e)] = n;
        }
        i = j + 1;
    }
    return m;
}

// Рамки стыковки узла: у сектора — криволинейные (грань ⊥ касательной опорной дуги), иначе ⊥ оси ящика.
// Общая для buildLine и BSCALE (смена угла полюса на ходу — рамки следуют за опорной).
inline void setFrames(Node& nd) {
    double Leff = nd.phys.Leff;
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
}

// ---- СБОРКА КОЛЬЦА из БД (build_line + check_junctions) ----------------------------------------
// nsl/order — срезы и порядок сектор/кик-тела (по умолчанию как elegant: 5 и 4); isr — квантовое возбуждение.
inline std::vector<Node> buildLine(const std::string& ws, int lattice, double design_mev, bool sr,
                                   const Settings& stg, int nsl_sector = SECTOR_NSLICES, int yoshida_order = 4,
                                   bool isr = false) {
    const std::map<std::string, int> nsl_dev = nslByDev();
    Table lc = readTable(ws + "tables/lattice_configs.csv");
    Table cn = readTable(ws + "tables/element_connectivity.csv");
    int cLat = lc.col("lattice_id"), cInst = lc.col("instance_id");
    int cx = lc.col("x"), cy = lc.col("y"), cz = lc.col("z"), cyaw = lc.col("yaw"),
        cpit = lc.col("pitch"), crol = lc.col("roll"), cdev = lc.col("device_type_id"),
        cB = lc.col("mat_B"), cG = lc.col("mat_G"), cGv = lc.col("mat_Gv"), cS = lc.col("mat_S"),
        cO = lc.col("mat_O"), cLeff = lc.col("eff_length"), cKind = lc.col("kind"),
        cNsl = lc.col("n_slices"), cE1 = lc.col("edge1"), cE2 = lc.col("edge2"),
        cLf1 = lc.col("lf1"), cLf2 = lc.col("lf2");             // L_f торцов (фриндж); нет столбца ⇒ 0
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
        { auto it = nsl_dev.find(nd.dev); if (it != nsl_dev.end()) nsl = it->second; }   // SS_NSL_DEV по типу
        double e1 = cellD(*r0, cE1), e2 = cellD(*r0, cE2);
        double lf1 = cellD(*r0, cLf1), lf2 = cellD(*r0, cLf2);
        int nsl_box = (int)cellD(*r0, cNsl, 40.0);
        nd.phys = Element::make(nd.kind, Leff, B, G, Gv, S, O, e1, e2, nsl, design_mev, sr, yoshida_order, isr,
                                nsl_box > 0 ? nsl_box : 40, lf1, lf2);
        nd.sector = nd.phys.sector; nd.monitor = (nd.kind == "monitor");
        double yaw = cellD(*r0, cyaw), pitch = cellD(*r0, cpit), roll = cellD(*r0, crol);
        nd.center = { cellD(*r0, cx), cellD(*r0, cy), cellD(*r0, cz) };
        frame(yaw, pitch, roll, nd.axis, nd.uhat, nd.vhat);
        nd.entry = nd.center - nd.axis * (0.5 * Leff);
        nd.exit  = nd.center + nd.axis * (0.5 * Leff);
        nd.theta0_0 = nd.phys.theta0; nd.edge1_0 = nd.phys.edge1; nd.edge2_0 = nd.phys.edge2;
        setFrames(nd);
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

// ---- ОПЫТ: пучок в рамке старт-плоскости, явные плоскости, резонаторы установки ----------------
// Порт beam_injection (Python-эталон): particle_groups — локальные (x, xp, y, yp, ct, delta) в рамке старта,
// experiments — start_x..start_nz / finish_x..finish_nz (пустой финиш ⇒ кольцо), rf_cavities — плоскости ВЧ.
struct LocalPart { double x = 0, xp = 0, y = 0, yp = 0, ct = 0, delta = 0, m = 0, q = 0; };

inline std::vector<LocalPart> loadGroup(const std::string& ws, int group) {
    Table pg = readTable(ws + "tables/particle_groups.csv");
    if (pg.col("xp") < 0)
        throw std::runtime_error("particle_groups.csv в старом формате — запустите migrate_explicit_planes.py");
    int g = pg.col("group_id"), cx = pg.col("x"), cxp = pg.col("xp"), cy = pg.col("y"), cyp = pg.col("yp"),
        cct = pg.col("ct"), cd = pg.col("delta"), cm = pg.col("m"), cq = pg.col("q");
    std::vector<LocalPart> out;
    for (auto& r : pg.rows) if (cell(r, g) == std::to_string(group)) {
        LocalPart p;
        p.x = cellD(r, cx); p.xp = cellD(r, cxp); p.y = cellD(r, cy); p.yp = cellD(r, cyp);
        p.ct = cellD(r, cct); p.delta = cellD(r, cd); p.m = cellD(r, cm); p.q = cellD(r, cq);
        out.push_back(p);
    }
    return out;
}

// (x, xp, y, yp, ct, δ) → канонические (x, px, y, py, t, δ): px = (1+δ)·xp/√(1+xp²+yp²) — local_to_canonical.
inline State localToCanonical(const LocalPart& p) {
    double k = (1.0 + p.delta) / std::sqrt(1.0 + p.xp * p.xp + p.yp * p.yp);
    return { p.x, k * p.xp, p.y, k * p.yp, p.ct, p.delta };
}

inline double gamma_of(double v) { return 1.0 / std::sqrt(1.0 - (v / C_LIGHT) * (v / C_LIGHT)); }
// опорный импульс P0 = bm.P_MOM = γ₀·m_e·v₀ (сурогат) — как bip_matrix.
inline double p0_ref() { return gamma_of(V_BEAM) * ME * V_BEAM; }

struct PlaneSpec { Vec3 c, n; };

// Плоскость из 6 ячеек (точка, нормаль): все пусты ⇒ false; заполнены частично / нулевая нормаль ⇒ ошибка.
inline bool readPlane(const Table& t, const std::vector<std::string>& r, const std::string& pre,
                      const char* const cols[6], PlaneSpec& out) {
    double v[6]; int filled = 0;
    for (int k = 0; k < 6; ++k) {
        std::string s = cell(r, t.col(pre + cols[k]));
        if (!s.empty() && s != "-") { v[k] = std::stod(s); ++filled; } else v[k] = 0.0;
    }
    if (filled == 0) return false;
    if (filled != 6) throw std::runtime_error("плоскость " + pre + ": заполнены не все 6 столбцов");
    out.c = { v[0], v[1], v[2] }; out.n = { v[3], v[4], v[5] };
    if (!(norm(out.n) > 0)) throw std::runtime_error("плоскость " + pre + ": нулевая нормаль");
    return true;
}

struct ExpPlanes { PlaneSpec start, finish; bool has_finish = false; };

inline ExpPlanes readExpPlanes(const Table& ex, const std::vector<std::string>& r) {
    static const char* const C6[6] = { "x", "y", "z", "nx", "ny", "nz" };
    if (ex.col("start_x") < 0)
        throw std::runtime_error("experiments.csv в старом формате (p1..p3/dir) — запустите migrate_explicit_planes.py");
    ExpPlanes p;
    if (!readPlane(ex, r, "start_", C6, p.start))
        throw std::runtime_error("опыт: не задана старт-плоскость (start_x..start_nz) — exp_planes.py --start-s 0");
    p.has_finish = readPlane(ex, r, "finish_", C6, p.finish);
    return p;
}

struct CavitySpec { std::string id; PlaneSpec pl; double share = 1.0; };

// Резонаторы установки (rf_cavities.csv; нет файла ⇒ нет резонаторов). Пустой volt_share ⇒ 1/N.
inline std::vector<CavitySpec> readCavities(const std::string& ws, int lattice) {
    std::vector<CavitySpec> out;
    std::ifstream probe(ws + "tables/rf_cavities.csv");
    if (!probe) return out;
    probe.close();
    static const char* const C6[6] = { "x", "y", "z", "nx", "ny", "nz" };
    Table t = readTable(ws + "tables/rf_cavities.csv");
    int ci = t.col("cavity_id"), cl = t.col("lattice_id"), cs = t.col("volt_share");
    std::vector<const std::vector<std::string>*> mine;
    for (auto& r : t.rows) if (cell(r, cl) == std::to_string(lattice)) mine.push_back(&r);
    for (auto* r : mine) {
        CavitySpec c; c.id = cell(*r, ci);
        if (!readPlane(t, *r, "", C6, c.pl)) throw std::runtime_error("rf_cavities: у резонатора " + c.id + " нет плоскости");
        std::string s = cell(*r, cs);
        c.share = (s.empty() || s == "-") ? 1.0 / (double)mine.size() : std::stod(s);
        out.push_back(c);
    }
    return out;
}

// Положение плоскости в линии — beam_injection.locate_plane: slot k — промежуток ПОСЛЕ mags[k] (N−1 — шов
// или луч после последнего), −1 — луч ДО первого; along — путь от выхода mags[k] (у −1 — от входа mags[0]).
inline void locatePlane(const std::vector<Node>& mags, const PlaneSpec& pl, int& slot, double& along) {
    Vec3 n = normalize(pl.n), c = pl.c;
    int N = (int)mags.size(); bool found = false; double bdist = 0;
    for (int i = 0; i < N; ++i) {
        Vec3 A = mags[i].exit, B = mags[(i + 1) % N].entry, AB = B - A;
        double L = norm(AB);
        if (L < 1e-9) continue;
        Vec3 axis = AB * (1.0 / L);
        double denom = dot(n, AB);
        if (std::fabs(denom) < 1e-12 || dot(n, axis) <= 0) continue;
        double t = dot(n, c - A) / denom;
        if (t >= -1e-9 && t <= 1.0 + 1e-9) {
            double dist = norm(A + AB * t - c);
            if (!found || dist < bdist) { found = true; bdist = dist; slot = i; along = t * L; }
        }
    }
    if (!found) {                                                  // открытая линия: лучи
        struct Ray { int slot; Vec3 P, ax; int sign; };
        Ray rays[2] = { { -1, mags[0].entry, mags[0].axis, -1 }, { N - 1, mags[N - 1].exit, mags[N - 1].axis, +1 } };
        for (auto& ry : rays) {
            double d = dot(n, ry.ax);
            if (d <= 1e-12) continue;
            double s = dot(n, c - ry.P) / d;
            if (ry.sign * s >= 0) {
                double dist = norm(ry.P + ry.ax * s - c);
                if (!found || dist < bdist) { found = true; bdist = dist; slot = ry.slot; along = s; }
            }
        }
    }
    if (!found) throw std::runtime_error("плоскость не пересекает опорную орбиту по направлению нормали "
                                         "ни в одном промежутке (внутри элемента или против хода?)");
}

// Линия — кольцо: шов (выход последнего → вход первого) — прямой промежуток вдоль их осей (beam_injection.SEAM_TOL).
inline bool seamIsDrift(const std::vector<Node>& mags) {
    const double TOL = 1e-3;
    Vec3 A = mags.back().exit, B = mags.front().entry, AB = B - A;
    double L = norm(AB);
    if (L < 1e-9) return true;
    Vec3 ax = AB * (1.0 / L);
    return dot(ax, mags.back().axis) > 0 && norm(cross(ax, mags.back().axis)) < TOL
        && norm(cross(ax, mags.front().axis)) < TOL;
}

// Совпадение плоскостей — matrix_tracker._plane_relation(...) == 'same' и нормали сонаправлены.
inline bool samePlane(const PlaneSpec& p, const PlaneSpec& q, const Settings& stg) {
    Vec3 nA = normalize(p.n), nB = normalize(q.n);
    double gap = std::fabs(dot(nB, p.c - q.c));
    return norm(cross(nA, nB)) < stg.plane_parallel_tol && gap < stg.plane_same_tol_m && dot(nA, nB) > 0;
}

// Результат сборки опыта: старт, путь (элементы + узлы-резонаторы), финиш, канон. состояния, кольцо/линия,
// опорная длина пролёта, дуги входов. end_share — доля резонаторов в плоскости финиша (< 0 ⇒ нет).
struct Injection { Node plane; std::vector<Node> ordered; Node finish; std::vector<State> states;
                   bool closed; double ring_len = 0.0; std::vector<double> entry_arc;
                   double end_share = -1.0; int n_cav = 0; std::vector<std::string> skipped; };

// Путь опыта без пучка — beam_injection.build_path: старт/финиш, элементы (с переходом через шов кольца),
// резонаторы cav в своих промежутках или в плоскости финиша. rf ⇒ при пустом cav ВЧ в конце пролёта.
inline Injection buildPath(const std::vector<Node>& mags, const ExpPlanes& ep, const std::vector<CavitySpec>& cav,
                           const Settings& stg, bool rf) {
    if (mags.empty()) throw std::runtime_error("решётка пуста");
    int N = (int)mags.size();
    Injection inj;
    inj.plane = makePlane(ep.start.c, ep.start.n);
    int ks = 0; double as = 0; locatePlane(mags, ep.start, ks, as);
    inj.closed = !ep.has_finish || samePlane(ep.start, ep.finish, stg);
    std::vector<int> order; int kf = ks; double af = as;
    if (inj.closed) {
        if (!seamIsDrift(mags)) throw std::runtime_error("финиш не задан (кольцо), но линия решётки не замкнута");
        inj.finish = inj.plane;
        for (int i = 0; i < N; ++i) order.push_back(((ks + 1 + i) % N + N) % N);
    } else {
        inj.finish = makePlane(ep.finish.c, ep.finish.n);
        locatePlane(mags, ep.finish, kf, af);
        if (kf > ks || (kf == ks && af > as)) {
            for (int i = ks + 1; i <= kf; ++i) order.push_back(i);
        } else {
            if (!seamIsDrift(mags)) throw std::runtime_error("финиш раньше старта, а линия не замкнута");
            for (int i = ks + 1; i < N; ++i) order.push_back(i);
            for (int i = 0; i <= kf; ++i) order.push_back(i);
        }
    }
    if (order.empty()) throw std::runtime_error("на пути старт → финиш нет ни одного элемента");
    // резонаторы установки (только с ВЧ): в своих промежутках пути или в плоскости финиша (конец пролёта)
    std::vector<int> slot_of_gap; slot_of_gap.push_back(ks);
    for (int i : order) slot_of_gap.push_back(i);
    size_t G = order.size();
    std::vector<std::vector<std::pair<double, Node>>> ins(G + 1);
    if (rf && cav.empty()) inj.end_share = 1.0;                        // прежнее правило: ВЧ в конце пролёта
    for (auto& cv : cav) {
        if (samePlane(cv.pl, inj.closed ? ep.start : ep.finish, stg)) {
            inj.end_share = (inj.end_share < 0 ? 0.0 : inj.end_share) + cv.share; continue;
        }
        if (!inj.closed && samePlane(cv.pl, ep.start, stg)) { inj.skipped.push_back(cv.id); continue; }
        int kc = 0; double ac = 0; locatePlane(mags, cv.pl, kc, ac);
        long hit = -1;
        for (size_t g = 0; g <= G; ++g) {
            if (slot_of_gap[g] != kc) continue;
            double lo = (g == 0) ? as : -INFINITY, hi = (g == G) ? af : INFINITY;
            if (lo < ac && ac < hi) { hit = (long)g; break; }
        }
        if (hit < 0) { inj.skipped.push_back(cv.id); continue; }
        Node nd = makePlane(cv.pl.c, cv.pl.n);
        nd.cavity = true; nd.cav_share = cv.share; nd.dev = "RF_CAVITY_" + cv.id;
        ins[hit].push_back({ ac, nd });
    }
    if (rf && !cav.empty() && inj.end_share < 0) {
        bool any = false; for (auto& v : ins) if (!v.empty()) any = true;
        if (!any) throw std::runtime_error("ни один резонатор установки (rf_cavities.csv) не стоит на пути старт → финиш");
    }
    for (size_t g = 0; g <= G; ++g) {
        std::stable_sort(ins[g].begin(), ins[g].end(),
                         [](const std::pair<double, Node>& a, const std::pair<double, Node>& b) { return a.first < b.first; });
        for (auto& pr : ins[g]) { inj.ordered.push_back(pr.second); ++inj.n_cav; }
        if (g < G) inj.ordered.push_back(mags[order[g]]);
    }
    // опорная арк-длина входа каждого узла (entry_arc) + длина пролёта (ring_len) — как beam_injection.build_path.
    auto& O = inj.ordered;
    inj.entry_arc.resize(O.size());
    inj.entry_arc[0] = dot(O[0].fr_in.p - inj.plane.fr_out.p, inj.plane.fr_out.axis);   // gap_in
    if (O[0].cavity) O[0].s_c = inj.entry_arc[0];
    for (size_t j = 1; j < O.size(); ++j) {
        double ord_gap = dot(O[j].fr_in.p - O[j - 1].fr_out.p, O[j - 1].fr_out.axis);
        inj.entry_arc[j] = inj.entry_arc[j - 1] + O[j - 1].phys.Leff + ord_gap;
        if (O[j].cavity) O[j].s_c = inj.entry_arc[j];
    }
    double gap_out = dot(inj.finish.fr_in.p - O.back().fr_out.p, O.back().fr_out.axis);
    inj.ring_len = inj.entry_arc.back() + O.back().phys.Leff + gap_out;
    return inj;
}

// Опыт → путь + канонические состояния пучка (rf ⇒ резонаторы установки встают в путь).
inline Injection inject(const std::string& ws, int lattice, int group, const ExpPlanes& ep, double design_mev,
                        bool sr, const Settings& stg, bool rf = false) {
    std::vector<Node> mags = buildLine(ws, lattice, design_mev, sr, stg);
    std::vector<LocalPart> parts = loadGroup(ws, group);
    if (parts.empty()) throw std::runtime_error("beam group not found");
    Injection inj = buildPath(mags, ep, rf ? readCavities(ws, lattice) : std::vector<CavitySpec>(), stg, rf);
    for (auto& p : parts) inj.states.push_back(localToCanonical(p));
    return inj;
}

// ---- ВЧ опыта: none | ideal | sine; резонаторы — узлы пути (cavity, s_c, cav_share) и/или конец пролёта ----
// Часы (дословный порт matrix_traj_export.SineRF): внутри пролёта t — «сырой» путь; на резонаторе запаздывание
// τ = t0 + Σ Δt_raw/β(δ) − (s_c/L)·C_RF, ΔE = V·share·sin(φ + 2πh·τ/C_RF); в конце пролёта
// t ← t0 + Σ Δt_raw/β(δ) − C_RF. Без резонаторов внутри — прежнее t ← t0 + путь/β(δ) − C_RF.
// ideal: δ ← ideal_delta — заданная энергия после резонатора (radiation_params.rf_ideal_energy_mev; пусто — проектная,
// δ = 0), какой бы частица ни прилетела (на каждом резонаторе и в конце, если там резонатор).
struct RfClock { double t0 = 0, raw = 0, acc = 0; bool fired = false;
                 void start(double t) { t0 = t; raw = t; acc = 0.0; fired = false; } };

struct RfSpec {
    enum Type { NONE = 0, IDEAL = 1, SINE = 2 } type = NONE;
    double volt = 0, h = 0, design = 6000.0, c_rf = 0, phase = 0, U0 = 0, L = 0, end_share = -1.0;
    double ideal_delta = 0.0;                        // ideal: δ после резонатора = E_ideal/E_design − 1
    bool auto_phase = false;
    const char* name() const { return type == IDEAL ? "ideal" : (type == SINE ? "sine" : "none"); }
    double beta(double d) const { double pc = design * (1.0 + d); return pc / std::sqrt(pc * pc + MC2_MEV * MC2_MEV); }
    void prepare(const Injection& inj, double freq_hz, bool phase_auto, double phase_deg);
    double kick(double d, double tau, double share) const {
        const double PI = 3.14159265358979323846;
        double dE = (volt * share) * std::sin(phase + 2.0 * PI * h * tau / c_rf);
        double pc = design * (1.0 + d);
        double E = std::sqrt(pc * pc + MC2_MEV * MC2_MEV) + dE;
        return std::sqrt(E * E - MC2_MEV * MC2_MEV) / design - 1.0;
    }
    void atCavity(State& s, const Node& c, RfClock& k) const {
        if (type == IDEAL) { s.delta = ideal_delta; return; }
        if (type != SINE) return;
        k.acc = k.acc + (s.t - k.raw) / beta(s.delta); k.raw = s.t; k.fired = true;
        double tau = k.t0 + k.acc - (c.s_c / L) * c_rf;
        s.delta = kick(s.delta, tau, c.cav_share);
    }
    void endPass(State& s, RfClock& k) const {
        if (type == IDEAL) { if (end_share > 0) s.delta = ideal_delta; return; }
        if (type != SINE) return;
        if (k.fired) s.t = k.t0 + (k.acc + (s.t - k.raw) / beta(s.delta)) - c_rf;
        else s.t = k.t0 + (s.t - k.t0) / beta(s.delta) - c_rf;
        if (end_share > 0) s.delta = kick(s.delta, s.t, end_share);
    }
};

// Один пролёт старт → путь → финиш; rf != nullptr ⇒ кики резонаторов на пути (конец пролёта — RfSpec::endPass).
inline State lineTurn(const Injection& inj, State s, const RfSpec* rf = nullptr, RfClock* k = nullptr) {
    const Frame* prev = &inj.plane.fr_out;
    for (const Node& nd : inj.ordered) {
        s = modified_drift(s, *prev, nd.fr_in, nd.same_plane_as_prev);
        s = nd.body(s);
        if (nd.cavity && rf) rf->atCavity(s, nd, *k);
        prev = &nd.fr_out;
    }
    return modified_drift(s, *prev, inj.finish.fr_in, inj.finish.same_plane_as_prev);
}

// Длина ВЧ-оборота и (при auto) синхронная фаза φs = π − asin(U0/V) по потере опорной частицы (нули).
inline void RfSpec::prepare(const Injection& inj, double freq_hz, bool phase_auto, double phase_deg) {
    const double PI = 3.14159265358979323846;
    L = inj.ring_len; end_share = inj.end_share;
    c_rf = (freq_hz > 0) ? h * C_LIGHT / freq_hz : inj.ring_len;
    auto_phase = phase_auto;
    if (phase_auto) {
        State z{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        double d1 = lineTurn(inj, z).delta;
        double pc0 = design, pc1 = design * (1.0 + d1);
        U0 = std::sqrt(pc0 * pc0 + MC2_MEV * MC2_MEV) - std::sqrt(pc1 * pc1 + MC2_MEV * MC2_MEV);
        if (!(std::fabs(U0) < volt)) throw std::runtime_error("ВЧ sine: напряжения не хватает на потерю за оборот");
        phase = PI - std::asin(U0 / volt);
    } else {
        phase = phase_deg * PI / 180.0;
    }
}

// ---- ТРЕК: N пролётов (старт → путь → финиш; у разомкнутой линии — ПОРТАЛ: финиш → старт) + ВЧ -------
struct TrackResult { State final_state; int surv_turns; bool lost; };

inline TrackResult trackParticle(const Injection& inj, const State& st0, int n_turns, const RfSpec& rf) {
    State s = st0; int done = 0; RfClock k;
    try {
        for (int turn = 0; turn < n_turns; ++turn) {
            k.start(s.t);
            s = lineTurn(inj, s, &rf, &k);
            rf.endPass(s, k);
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

inline LongResult trackLong(const Injection& inj, const State& st0, int n_turns, const RfSpec& rf) {
    LongResult r; State s = st0; RfClock k;
    r.tt.push_back({ s.x, s.px, s.y, s.py });
    try {
        for (int turn = 0; turn < n_turns; ++turn) {
            k.start(s.t);
            s = lineTurn(inj, s, &rf, &k);
            rf.endPass(s, k);
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
            // финиш — явная плоскость опыта (кольцо: = старт, элемент 0; линия: элемент N+1, дальше портал)
            st = modified_drift(st, O.back().fr_out, inj.finish.fr_in, inj.finish.same_plane_as_prev);
            rec(base + inj.ring_len, st, inj.closed ? 0 : (int64_t)(O.size() + 1), inj.finish.fr_in);
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
