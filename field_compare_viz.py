# -*- coding: utf-8 -*-
"""Наглядное сравнение РЕАЛЬНОГО поля (карта ANSYS, его видит Boris) и АНАЛИТИЧЕСКОГО
идеального мультиполя (его используют elegant/MAD/BDSIM) — для ВСЕХ элементов.
Идеал = чистый мультиполь порядка N с комплексной силой, подогнанной к карте; разница
(реал − идеал) = «неидеальность» магнита (высшие мультиполи, краевые поля).

Фигуры (workspace/bip_compare/figs/):
  fig_realvideal_<dev>.png — поперечное поле (реал / идеал / разница ×k), тепловые карты;
  fig_field_radial.png     — сетка: |B⊥| вдоль радиуса, реал vs идеал ∝ r^(N-1);
  fig_field_fringe.png     — сетка: продольный профиль силы (fringe vs hard-edge).
"""
import struct, os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HDR = 96
FIELD_DIR = os.path.join("workspace", "fields_bin")
FIGS = os.path.join("workspace", "bip_compare", "figs")
# (device, заголовок, orders — порядки мультиполя в «идеальной» модели, main — главный для срезов).
# Чистые элементы: один порядок. Диполь-квадруполи combined-function: диполь+квадруполь (1,2).
ELEMS = [
    ("QP_L_1A", "Квадруполь длинный", (2,), 2),
    ("QP_sh_1A", "Квадруполь короткий", (2,), 2),
    ("SP_I1_1A", "Секступоль (гл. обм.)", (3,), 3),
    ("SP_I2_1A", "Секступоль (корр. обм.)", (3,), 3),
    ("OP_1A", "Октуполь", (4,), 4),
    ("DQ_I1_1A", "Диполь-квадруполь I1", (1, 2), 1),
    ("DQ_I2_1A", "Диполь-квадруполь I2", (1, 2), 2),
    ("DQ_I3_1A", "Диполь-квадруполь I3", (1, 2), 2),
    ("DP_L_PM", "Диполь длинный", (1,), 1),
    ("DP_sh_PM", "Диполь короткий", (1,), 1),
    ("corrector_I1_1A", "Корректор I1", (1,), 1),
    ("corrector_I2_1A", "Корректор I2", (1,), 1),
]


def load(device):
    with open(os.path.join(FIELD_DIR, device + ".bin"), "rb") as f:
        blob = f.read(HDR)
        mnx, mny, mnz, mxx, mxy, mxz, sx, sy, sz = struct.unpack_from("<9d", blob, 0)
        nx, ny, nz = struct.unpack_from("<3i", blob, 72)
        B = np.frombuffer(f.read(nx * ny * nz * 24), dtype="<f8").reshape(nx, ny, nz, 3)
    return dict(mny=mny, mnz=mnz, sx=sx, sy=sy, sz=sz, nx=nx, ny=ny, nz=nz), B


def peak_slice(h, B):
    j0, k0 = h["ny"] // 2, h["nz"] // 2
    mag = np.sqrt(np.sum(B[:, j0, k0] ** 2, axis=-1))
    return int(np.argmax(mag)) if mag.max() > 0 else 0


def fit_coeffs(Y, Z, By, Bz, orders, mask):
    """МНК-подгон набора мультиполей: B_z + iB_y = Σ K_n·(y+iz)^(n-1) (K комплексные)."""
    cols = [((Y + 1j * Z) ** (n - 1))[mask].ravel() for n in orders]
    A = np.column_stack(cols)
    b = (Bz + 1j * By)[mask].ravel()
    K, *_ = np.linalg.lstsq(A, b, rcond=None)
    return dict(zip(orders, K))


def ideal_field(Y, Z, coeffs):
    f = sum(K * (Y + 1j * Z) ** (n - 1) for n, K in coeffs.items())
    return f.imag, f.real          # By, Bz


def fig_transverse(dev, title, orders, n):
    h, B = load(dev)
    ix = peak_slice(h, B)
    j0, k0 = h["ny"] // 2, h["nz"] // 2
    R = 8
    js = slice(j0 - R, j0 + R + 1); ks = slice(k0 - R, k0 + R + 1)
    y = (np.arange(h["ny"])[js] - j0) * h["sy"]; z = (np.arange(h["nz"])[ks] - k0) * h["sz"]
    Z, Ym = np.meshgrid(z, y)
    By = B[ix, js, ks, 1]; Bz = B[ix, js, ks, 2]
    rr = np.sqrt(Ym ** 2 + Z ** 2); mask = rr <= 5e-3
    coeffs = fit_coeffs(Ym, Z, By, Bz, orders, mask)
    iBy, iBz = ideal_field(Ym, Z, coeffs)
    # доминирующая компонента (для диполей — равномерное поле в своей оси)
    if np.abs(Bz).max() >= np.abs(By).max():
        real, ideal, clab = Bz, iBz, "Bz"
    else:
        real, ideal, clab = By, iBy, "By"
    res = real - ideal
    ext = [z[0] * 1e3, z[-1] * 1e3, y[0] * 1e3, y[-1] * 1e3]
    vmax = np.abs(real).max() * 1e6 + 1e-9
    rmax = np.abs(res).max() * 1e6 + 1e-9
    kamp = max(1, int(round(vmax / rmax / 1.5)))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.7))
    data = [(f"Реальное поле {clab} (ANSYS)", real * 1e6, vmax),
            (f"Идеальный мультиполь {clab}", ideal * 1e6, vmax),
            (f"Разница (реал − идеал) ×{kamp}", res * 1e6 * kamp, vmax)]
    for ax, (ttl, D, vm) in zip(axes, data):
        im = ax.imshow(D, origin="lower", extent=ext, cmap="RdBu_r", vmin=-vm, vmax=vm, aspect="equal")
        ax.add_patch(plt.Circle((0, 0), 2, fill=False, ls="--", color="k", lw=1.2))
        ax.set_xlabel("z, мм"); ax.set_ylabel("y, мм"); ax.set_title(ttl, fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=f"{clab}, мкТл")
    nonideal = math.sqrt(np.mean((By - iBy)[mask] ** 2 + (Bz - iBz)[mask] ** 2)) / \
        (math.sqrt(np.mean(By[mask] ** 2 + Bz[mask] ** 2)) + 1e-30)
    omod = "+".join(str(o) for o in orders)
    fig.suptitle(f"{title} ({dev}): поперечное поле в центре магнита (идеал = мультип. N={omod}). "
                 f"Пунктир — пучок r=2 мм. Отклонение от идеала (r≤5мм) ≈ {nonideal*100:.0f}%", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = os.path.join(FIGS, f"fig_realvideal_{dev}.png")
    fig.savefig(p, dpi=110); plt.close(fig); return p, nonideal


def _bperp_ray(h, B, ix, r, th):
    y, z = r * math.cos(th), r * math.sin(th)
    j0, k0 = h["ny"] // 2, h["nz"] // 2
    fy = y / h["sy"] + j0; fz = z / h["sz"] + k0
    j, k = int(fy), int(fz); ty, tz = fy - j, fz - k
    out = []
    for c in (1, 2):
        g = B[ix, :, :, c]
        out.append(g[j, k] * (1 - ty) * (1 - tz) + g[j + 1, k] * ty * (1 - tz) +
                   g[j, k + 1] * (1 - ty) * tz + g[j + 1, k + 1] * ty * tz)
    return math.hypot(out[0], out[1])


def fig_radial_grid():
    ncol = 4; nrow = math.ceil(len(ELEMS) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3 * nrow))
    axes = np.array(axes).reshape(-1)
    for ax, (dev, title, orders, n) in zip(axes, ELEMS):
        h, B = load(dev); ix = peak_slice(h, B)
        j0, k0 = h["ny"] // 2, h["nz"] // 2; R = 6
        js = slice(j0 - R, j0 + R + 1); ks = slice(k0 - R, k0 + R + 1)
        yg = (np.arange(h["ny"])[js] - j0) * h["sy"]; zg = (np.arange(h["nz"])[ks] - k0) * h["sz"]
        Zg, Yg = np.meshgrid(zg, yg)
        rrg = np.sqrt(Yg ** 2 + Zg ** 2)
        coeffs = fit_coeffs(Yg, Zg, B[ix, js, ks, 1], B[ix, js, ks, 2], orders, rrg <= 5e-3)
        th = math.radians(30); rs = np.linspace(0, 8e-3, 40)
        real = np.array([_bperp_ray(h, B, ix, r, th) for r in rs])
        fid = sum(K * (rs * math.cos(th) + 1j * rs * math.sin(th)) ** (o - 1) for o, K in coeffs.items())
        ideal = np.abs(fid)
        olab = "∝r^%d" % (n - 1) if len(orders) == 1 else "combined " + "+".join(map(str, orders))
        ax.plot(rs * 1e3, real * 1e6, "o-", ms=2.5, color="tab:blue", label="реал")
        ax.plot(rs * 1e3, ideal * 1e6, "--", color="tab:green", label=f"идеал {olab}")
        ax.axvspan(0, 2, color="0.9")
        ax.set_title(title, fontsize=8.5); ax.grid(True, alpha=0.3); ax.tick_params(labelsize=7)
        ax.set_xlabel("r, мм", fontsize=8); ax.set_ylabel("|B⊥|, мкТл", fontsize=8); ax.legend(fontsize=6.5)
    for ax in axes[len(ELEMS):]:
        ax.axis("off")
    fig.suptitle("Поле вдоль радиуса (под 30°): реальное (синий) vs идеальный мультиполь (зелёный). "
                 "Серое — зона пучка. Отклонение = высшие мультиполи / краевые поля", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    p = os.path.join(FIGS, "fig_field_radial.png"); fig.savefig(p, dpi=105); plt.close(fig); return p


def fig_fringe_grid():
    ncol = 4; nrow = math.ceil(len(ELEMS) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 2.7 * nrow))
    axes = np.array(axes).reshape(-1)
    for ax, (dev, title, orders, n) in zip(axes, ELEMS):
        h, B = load(dev)
        j0, k0 = h["ny"] // 2, h["nz"] // 2; R = 5
        js = slice(j0 - R, j0 + R + 1); ks = slice(k0 - R, k0 + R + 1)
        y = (np.arange(h["ny"])[js] - j0) * h["sy"]; z = (np.arange(h["nz"])[ks] - k0) * h["sz"]
        Z, Y = np.meshgrid(z, y); w = (Y + 1j * Z) ** (n - 1)
        f = (B[:, js, ks, 2] + 1j * B[:, js, ks, 1])
        Kx = np.abs((f * np.conj(w)).sum(axis=(1, 2)) / (np.sum(np.abs(w) ** 2) + 1e-300))
        ic = int(np.argmax(Kx)); pk = Kx.max() + 1e-30
        xrel = (np.arange(h["nx"]) - ic) * h["sx"] * 1e3
        mm = xrel >= 0
        # эфф. полудлина (hard-edge) = площадь/пик /2
        L2 = float(np.sum(Kx) * h["sx"] / pk) * 1e3 / 2
        ax.plot(xrel[mm], (Kx / pk)[mm], lw=1.6, color="tab:blue", label="реал (fringe)")
        ax.plot([0, L2, L2], [1, 1, 0], "k--", lw=1.3, label="идеал hard-edge")
        ax.set_title(title, fontsize=8.5); ax.grid(True, alpha=0.3); ax.tick_params(labelsize=7)
        ax.set_xlabel("|X−центр|, мм", fontsize=8); ax.set_ylabel("|K| норм.", fontsize=8)
        ax.set_xlim(0, max(80, L2 * 1.3)); ax.set_ylim(-0.03, 1.08); ax.legend(fontsize=6.5)
    for ax in axes[len(ELEMS):]:
        ax.axis("off")
    fig.suptitle("Продольный профиль силы: реальное поле — плавный краевой спад (fringe), "
                 "идеал — резкий обрыв (hard-edge) той же площади ∫K·dx", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    p = os.path.join(FIGS, "fig_field_fringe.png"); fig.savefig(p, dpi=105); plt.close(fig); return p


def main():
    os.makedirs(FIGS, exist_ok=True)
    print("Сравнение реального и идеального поля по всем элементам:")
    res = []
    for dev, title, orders, n in ELEMS:
        try:
            p, ni = fig_transverse(dev, title, orders, n)
            res.append((dev, title, orders, ni))
            print(f"  {title:<30} N={'+'.join(map(str,orders)):<4} откл. от идеала(r≤5мм) ≈ {ni*100:5.1f}%")
        except FileNotFoundError:
            print(f"  {title:<30} — карта не найдена, пропуск")
    print("  сетка радиальных срезов: " + os.path.basename(fig_radial_grid()))
    print("  сетка продольных профилей: " + os.path.basename(fig_fringe_grid()))
    import csv
    with open(os.path.join("workspace", "bip_compare", "field_nonideality.csv"), "w",
              newline="", encoding="utf-8") as fcsv:
        wr = csv.writer(fcsv, delimiter=";")
        wr.writerow(["device", "name", "ideal_orders", "deviation_pct_r5mm"])
        for dev, title, orders, ni in res:
            wr.writerow([dev, title, "+".join(map(str, orders)), f"{ni*100:.2f}"])


if __name__ == "__main__":
    main()
