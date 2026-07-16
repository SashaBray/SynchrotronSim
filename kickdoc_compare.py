# -*- coding: utf-8 -*-
"""Сравнение траекторий ДИФФУР (Boris по картам поля) vs КИК-ЯДРО (matrix_tracker) на КАЖДОМ
элементе библиотеки — для технического документа. probe-пучок смещённых частиц через один элемент;
метрики совпадения; коэффициенты кик-представления; матрица связи; таблица ток→физкоэф.
Данные → workspace/kickdoc/data.json, фигуры → workspace/kickdoc/figs/. Запуск: py -3.14 kickdoc_compare.py"""
import os, json, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matrix_tracker as mt
import dq_sector_profile as dsp
from bip_matrix import QP

OUT = os.path.join("workspace", "research", "kickdoc"); FIGS = os.path.join(OUT, "figs")
os.makedirs(FIGS, exist_ok=True)
plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans", "figure.dpi": 130})
HW = {r["device_type_id"]: r for r in mt._read("hardware_library.csv")}
CM = {}                                  # coupling_matrix УДАЛЁН (токовый режим вырезан) — coupling-колонки пусты
EFF, PHYS, COUP, TRANS = mt.load_db()

# (dev, значение-аргумента, подпись). Для current — ток [А]; для coeff (диполи) — correct_coeff.
# проектно-масштабные аргументы (малые углы, линейный режим — где кик-модель и есть приближение поля)
ELEMS = [
    ("QP_L_1A", 2.0, "Квадруполь длинный"), ("QP_sh_1A", 3.0, "Квадруполь короткий"),
    ("SP_I1_1A", 10.0, "Секступоль I1"), ("SP_I2_1A", 10.0, "Секступоль I2"),
    ("OP_1A", 10.0, "Октуполь"),
    ("DQ_I1_1A", 5.0, "Диполь-квадруполь, обмотка I1"), ("DQ_I2_1A", 5.0, "Диполь-квадруполь, обмотка I2"),
    ("DQ_I3_1A", 5.0, "Диполь-квадруполь, обмотка I3"),
    ("DP_L_PM", 1.0, "Поворотный ПМ длинный"), ("DP_sh_PM", 1.0, "Поворотный ПМ короткий"),
    ("corrector_I1_1A", 20.0, "Корректор горизонтальный"), ("corrector_I2_1A", 20.0, "Корректор вертикальный"),
]
PROBES = [(0.0, 0.0), (2e-3, 0.0), (-2e-3, 0.0), (0.0, 2e-3), (0.0, -2e-3)]  # (u0,v0), м


def scaled_field(dev, val):
    """Поле элемента при аргументе val (ток или коэф.), с корректными множителями. → field_fn, half."""
    h, B, nom = dsp._load_full(dev, HW[dev])
    itype = CM.get(dev, {}).get("input_type", "current")
    fd = float(HW[dev].get("field_direction") or 1)
    hwcc = float(HW[dev].get("correct_coeff") or 1)
    scale = (val / nom) if itype == "current" else (val * hwcc * fd)
    return dsp._field_fn(h, B * scale), float(h["mxx"])


def boris(field, half, u0, v0, ds=2e-4):
    """Boris по полю: probe (u0,v0,u'=0,v'=0) вдоль +X через зону [-half,half]. → s,u,v,up,vp."""
    x = -half; u, v, up, vp = u0, v0, 0.0, 0.0
    n = int(round(2 * half / ds)); S, U, V = [0.0], [u], [v]
    for _ in range(n):
        ps = math.sqrt(max(1e-12, 1 - up * up - vp * vp)); x += ds / 2; u += up * ds / (2 * ps); v += vp * ds / (2 * ps)
        B = field(np.array([x, u, v])); up += -QP * B[2] * ds; vp += QP * B[1] * ds
        ps = math.sqrt(max(1e-12, 1 - up * up - vp * vp)); x += ds / 2; u += up * ds / (2 * ps); v += vp * ds / (2 * ps)
        S.append(x + half); U.append(u); V.append(v)
    return np.array(S), np.array(U), np.array(V), up, vp


def kick_probe(m, half, u0, v0):
    """Кик-ядро (matrix_tracker.Magnet) над зоной поля 2·half (для совпадения зоны с Boris);
    сила кика по Leff (реальная). Возврат s,u,v + финальные up,vp. Реплика Magnet.body по срезам."""
    N = m.nsl
    ds_g = 2.0 * half / N; ds_f = m.Leff / N
    th, k1, k1v, k2, k3, kE = m.theta0, m.k1, m.k1v, m.k2, m.k3, m.kEdge
    st = (u0, 0.0, v0, 0.0, 0.0, 0.0); s = 0.0; S, U, Vv = [0.0], [u0], [v0]
    for _ in range(N):
        st = mt._drift(st, ds_g / 2.0); s += ds_g / 2.0
        x, px, y, py, t, d = st
        px += th / N - k1 * x * ds_f - 0.5 * k2 * (x * x - y * y) * ds_f - (k3 / 6.0) * (x**3 - 3 * x * y * y) * ds_f
        py += k1v * y * ds_f + k2 * x * y * ds_f + (k3 / 6.0) * (3 * x * x * y - y**3) * ds_f - (kE / N) * y
        st = mt._drift((x, px, y, py, t, d), ds_g / 2.0); s += ds_g / 2.0
        S.append(s); U.append(st[0]); Vv.append(st[2])
    return np.array(S), np.array(U), np.array(Vv), st[1], st[3]


def make_magnet(dev, val):
    itype = CM.get(dev, {}).get("input_type", "current")
    inst = {"device_type_id": dev, "x": "0", "y": "0", "z": "0", "yaw": "0", "pitch": "0", "roll": "0",
            "instance_id": "1", "arg_val": str(val) if itype == "current" else "0",
            "correct_coeff": str(val) if itype == "coeff" else "1"}
    return mt.Magnet([inst], EFF, PHYS, COUP, TRANS, source="current")


def coef_table(dev):
    """Таблица ток→физкоэф. по матрице связи (5 значений тока)."""
    c = CM.get(dev, {}); itype = c.get("input_type", "current")
    perB, perG, perS, perO = (float(c.get(k) or 0) for k in ("per_B_Tl", "per_G_Tlm", "per_S_Tlm2", "per_O_Tlm3"))
    vals = [-100, -50, 0, 50, 100] if itype == "current" else [-2, -1, 0, 1, 2]
    rows = []
    for I in vals:
        rows.append(dict(I=I, B=perB * I, G=perG * I, S=perS * I, O=perO * I))
    return dict(itype=itype, per=dict(B=perB, G=perG, S=perS, O=perO), rows=rows)


def main():
    data = {"QP": float(QP), "elements": []}
    for dev, val, title in ELEMS:
        m = make_magnet(dev, val)
        field, half = scaled_field(dev, val)
        # корр.I2 (kick_z) установлен в кольце с roll 90° → его карта (гориз. изгиб) в кольце
        # отклоняет вертикально. Разворачиваем диффур-результат на 90° (u,v)->(-v,u), чтобы
        # сравнивать в той же плоскости, что и вертикальный кик модели.
        roll90 = "kick_z" in (TRANS[dev].get("active_coeffs") or "")
        # probe-пучок: diff + kick
        met_up, met_vp = [], []; trajs = []
        for u0, v0 in PROBES:
            Sb, Ub, Vb, upb, vpb = boris(field, half, u0, v0)
            if roll90:
                Ub, Vb = -Vb, Ub.copy(); upb, vpb = -vpb, upb
            Sk, Uk, Vk, upk, vpk = kick_probe(m, half, u0, v0)
            met_up.append(abs(upk - upb)); met_vp.append(abs(vpk - vpb))
            trajs.append(dict(u0=u0, v0=v0, Sb=Sb.tolist(), Ub=Ub.tolist(), Vb=Vb.tolist(),
                              Sk=Sk.tolist(), Uk=Uk.tolist(), Vk=Vk.tolist()))
        # фигура: u(s) и v(s), diff vs kick, для смещённых probe
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
        cols = ["#000", "#c0392b", "#2980b9", "#27ae60", "#e67e22"]
        for i, tr in enumerate(trajs):
            ax[0].plot(tr["Sb"], np.array(tr["Ub"]) * 1e3, color=cols[i], lw=2.4, alpha=.5)
            ax[0].plot(tr["Sk"], np.array(tr["Uk"]) * 1e3, color=cols[i], lw=1.0, ls="--")
            ax[1].plot(tr["Sb"], np.array(tr["Vb"]) * 1e3, color=cols[i], lw=2.4, alpha=.5)
            ax[1].plot(tr["Sk"], np.array(tr["Vk"]) * 1e3, color=cols[i], lw=1.0, ls="--")
        ax[0].set_xlabel("s, м"); ax[0].set_ylabel("u (гор.), мм"); ax[0].set_title(f"{dev}: гориз."); ax[0].grid(alpha=.3)
        ax[1].set_xlabel("s, м"); ax[1].set_ylabel("v (верт.), мм"); ax[1].set_title("верт. (— диффур, -- кик)"); ax[1].grid(alpha=.3)
        fig.tight_layout(); fig.savefig(os.path.join(FIGS, f"{dev}.png")); plt.close(fig)
        data["elements"].append(dict(
            dev=dev, title=title, val=val, kind=m.kind, nsl=m.nsl, Leff=m.Leff, Lphys=m.L,
            input_type=CM.get(dev, {}).get("input_type") or "current", active=TRANS[dev].get("active_coeffs"),
            coeffs=dict(theta0=m.theta0, k1=m.k1, k2=m.k2, k3=m.k3, kEdge=m.kEdge),
            metric=dict(max_dup_urad=max(met_up) * 1e6, max_dvp_urad=max(met_vp) * 1e6),
            coupling=coef_table(dev)))
        print(f"  {dev:16} kind={m.kind:7} θ0={m.theta0*1e3:+.3f}мрад k1={m.k1:+.3f} k2={m.k2:+.2f} "
              f"| Δu'max={max(met_up)*1e6:.2f}мкрад Δv'max={max(met_vp)*1e6:.2f}мкрад")
    json.dump(data, open(os.path.join(OUT, "data.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"[OK] {len(data['elements'])} элементов → {OUT}/data.json + figs/")


if __name__ == "__main__":
    main()
