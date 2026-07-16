# -*- coding: utf-8 -*-
"""Валидация таблицы element_coupling.csv end-to-end через МАТРИЧНЫЙ трекер (не через field-Boris).
Для набора эталонных настроек корректора: суперпозиция мультиполей из таблицы × токи → матричный
Magnet → дрейф+тело+дрейф (геометрия эталона) → сравнение финалов с corr_traj.pkl. Если RMS сравним
с corrector_model_ho (~0.07мм) + тонколинзовая поправка — таблица корректна (units не перепутаны)."""
import os, csv, pickle, math
import numpy as np
import matrix_tracker as mt
import corrector_experiment as ce

MULT = ["per_B_Tl", "per_G_Tlm", "per_S_Tlm2", "per_O_Tlm3"]
MATCOL = ["mat_B", "mat_G", "mat_S", "mat_O"]
WIND = ["main", "c1", "c2", "c3", "c4"]
LEFF, LMET = 0.1615, 0.108


def load_matrix():
    rows = {r["winding"]: r for r in csv.DictReader(
        open(os.path.join("workspace", "tables", "element_coupling.csv"), encoding="utf-8"), delimiter=";")
        if r["element_id"] == "corrector"}
    return {w: [float(rows[w][m]) for m in MULT] for w in WIND}      # per-1А мультиполи по обмоткам


def main():
    M = load_matrix()
    ref = pickle.load(open(os.path.join(ce.OUT, "corr_traj.pkl"), "rb")); meta = ref["meta"]
    hs, Bs = zip(*[ce.load_cached(m) for m in ce.MAPS]); h = hs[0]
    half = h["mxx"] - h["mnx"]
    XSTART = meta["XSTART"]
    drift_a = (-LMET / 2) - XSTART                                   # старт → вход коробки (центр элемента X=0)
    drift_b = (half + meta["drift2"]) - (LMET / 2)                   # выход коробки → финиш
    print(f"Геометрия: half={half*1e3:.0f}мм, drift_a={drift_a*1e3:.0f}мм, коробка L={LMET*1e3:.0f}мм, drift_b={drift_b*1e3:.0f}мм")

    dev = "CT"; eff = {dev: LEFF}; ph = {dev: LMET}
    tr = {dev: {"kind": "multipole", "model": "thin", "n_slices": 40}}

    def make_mag(cur):
        sup = [sum(cur[j] * M[WIND[j]][i] for j in range(5)) for i in range(8)]   # суперпозиция × токи (А)
        inst = [dict({"device_type_id": dev, "x": 0, "y": 0, "z": 0, "yaw": 0},
                     **{MATCOL[i]: sup[i] for i in range(8)})]
        return mt.Magnet(inst, eff, ph, {}, tr, source="physical")

    def track(sett):
        mag = make_mag(sett["currents"]); errs = []
        for gname, g in sett["groups"].items():
            for i in range(len(g["init"])):
                u, up, v, vp = g["init"][i]
                try:
                    st = mt._drift([u, up, v, vp, 0.0, 0.0], drift_a)
                    st = mag.body(st); st = mt._drift(st, drift_b)
                except mt.ParticleLost:
                    continue
                rp = g["final_pos"][i]
                errs.append(math.hypot(st[0] - rp[1], st[2] - rp[2]))
        return errs

    SUB = ["main_+100", "correct_1_+100", "correct_2_+100", "quad_norm", "quad_skew",
           "rand23", "rand01", "rand09", "rand05", "rand00"]
    SUB = [s for s in SUB if s in ref["settings"]]
    print(f"\n{'настройка':>16} {'RMS,мм':>9} {'макс,мм':>9}  токи")
    alle = []
    for lab in SUB:
        sett = ref["settings"][lab]; e = track(sett); alle += e
        if e:
            print(f"{lab:>16} {np.sqrt(np.mean(np.square(e)))*1e3:>9.3f} {np.max(e)*1e3:>9.2f}  {sett['currents']}")
    print(f"\nОбщий RMS по {len(SUB)} настройкам: {np.sqrt(np.mean(np.square(alle)))*1e3:.3f}мм")
    print("(corrector_model_ho full-field: 0.070мм; матричная тонколинзовая коробка чуть грубее)")


if __name__ == "__main__":
    main()
