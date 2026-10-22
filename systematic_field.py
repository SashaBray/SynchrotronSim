# -*- coding: utf-8 -*-
"""ПУНКТ B статики: СИСТЕМАТИКА ПОЛЯ ДИПОЛЕЙ (когерентная ошибка семейства) и парирование ВЧ.

Случайная ошибка поля → произвольная орбита (корректоры). КОГЕРЕНТНАЯ (все диполи ×(1+δ_f)) → это
СДВИГ ЭНЕРГИИ: изгиб ∝ B/p, поэтому поле +δ_f ≡ импульс −δ_f. Пучок хранится на ДИСПЕРСИОННОЙ орбите
η·δ, и это убирается ВЧ-ЧАСТОТОЙ (Δf/f = −α·Δp/p), а не корректорами.

Проверяем: (1) масштабирование поля диполей ≡ off-momentum орбита (эквивалентность); (2) компакция α и
трим Δf/f для компенсации; (3) хроматический сдвиг тюна ξ·δ_f. Запуск: py -3.14 systematic_field.py
"""
import os, sys, csv, tempfile, shutil
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, ring_analysis as ra
from bip_matrix import QP

inf = dict(fr.info_of(26)); inf["_id"] = 26
r = fr._ring(inf, nsl=8)
C = float(r.circumference)
co0 = ra.closed_orbit(r, dims=4); o0 = r.track_elements([co0])[0]
base_rms = float(np.sqrt(np.mean(o0[1:, 0] ** 2))) * 1e6
# компакция α = (1/C)·Σ_диполей η_i·θ_i  (θ_i = −mat_B·QP·L из mat_of)
tw = ra.twiss(r); eta = np.array(tw["etax"])
th = np.array([-float(getattr(m, "mat_B", 0.0) or 0.0) * QP * float(getattr(m, "Leff", 0.0) or 0.0)
               for m in r.mags])
n = min(len(eta), len(th))
alpha = float(np.sum(eta[:n] * th[:n]) / C)
# хроматизм
chrom = ra.chromaticity(r, deltas=np.linspace(-5e-3, 5e-3, 7), order=2)
xix, xiy = float(chrom["xi_x"][1]), float(chrom["xi_y"][1])
print(f"СИСТЕМАТИКА ПОЛЯ ДИПОЛЕЙ (реш.26). Периметр {C:.2f} м, α={alpha:.3e}, ξx={xix:.1f} ξy={xiy:.1f}")
print(f"базовая орбита идеала: {base_rms:.2f} мкм\n")

print(f"{'δ_f поле':>9} | off-momentum δp=−δ_f | трим ВЧ Δf/f | сдвиг тюна νx,νy")
for df in (1e-4, 5e-4, 1e-3, 2e-3):
    co = ra.closed_orbit(r, delta=-df, dims=4)
    if not np.all(np.isfinite(co)):
        print(f"{df*100:7.3f}% | орбита не замкнулась"); continue
    o = r.track_elements([co])[0]
    rms = float(np.sqrt(np.mean(o[1:, 0] ** 2))) * 1e6; mx = float(np.max(np.abs(o[1:, 0]))) * 1e6
    dff = alpha * df                                      # Δf/f для компенсации (= α·Δp/p)
    print(f"{df*100:7.3f}% | орбита rms={rms:6.1f} max={mx:6.0f} мкм | {dff:+.3e} | {xix*df:+.4f},{xiy*df:+.4f}")

# ВЕРИФИКАЦИЯ: реально масштабируем mat_B диполей ×(1+δ_f) и сравниваем с off-momentum
DF = 5e-4
root = tempfile.mkdtemp(prefix="sysfld_")
TB = "workspace/tables"
rows = list(csv.reader(open(os.path.join(TB, "lattice_configs.csv"), encoding="utf-8-sig"), delimiter=";"))
h = rows[0]; ci = {c: i for i, c in enumerate(h)}
os.makedirs(os.path.join(root, "tables"))
for rr in rows[1:]:
    if rr and rr[ci["lattice_id"]] == "26" and rr[ci["kind"]].strip() in ("dipole", "dq"):
        for col in ("mat_B", "mat_G", "mat_S", "mat_O"):
            if rr[ci[col]] not in ("", "0"):
                rr[ci[col]] = repr(float(rr[ci[col]]) * (1 + DF))
with open(os.path.join(root, "tables", "lattice_configs.csv"), "w", newline="", encoding="utf-8") as f:
    csv.writer(f, delimiter=";").writerows(rows)
for t in ("element_connectivity.csv", "particle_groups.csv"):
    shutil.copyfile(os.path.join(TB, t), os.path.join(root, "tables", t))
r2 = fr._ring(inf, ws=root, nsl=8)
co_sf = ra.closed_orbit(r2, dims=4); o_sf = r2.track_elements([co_sf])[0]
rms_sf = float(np.sqrt(np.mean(o_sf[1:, 0] ** 2))) * 1e6
# off-momentum для сравнения
co_om = ra.closed_orbit(r, delta=-DF, dims=4); o_om = r.track_elements([co_om])[0]
rms_om = float(np.sqrt(np.mean(o_om[1:, 0] ** 2))) * 1e6
corr = float(np.corrcoef(o_sf[1:n, 0], o_om[1:n, 0])[0, 1])
print(f"\nВЕРИФИКАЦИЯ (δ_f={DF*100:.2f}%): поле×(1+δ_f) орбита rms={rms_sf:.1f} мкм vs off-momentum rms={rms_om:.1f} мкм, corr орбит={corr:.4f}")
r2.close(); r.close(); shutil.rmtree(root, ignore_errors=True)
print("\n[готово] систематика = энергосдвиг → ВЧ-частота; corr≈1 подтверждает эквивалентность")
