# -*- coding: utf-8 -*-
"""Срез поля вдоль произвольной прямой P1→P2: три графика Bx, By, Bz по длине линии.
Реальное поле (карта ANSYS, 3D-трилинейная выборка) против аналитического идеального
мультиполя (подогнан в центре магнита). Точки — в локальных координатах карты (мм):
X — продольная ось, Y и Z — поперечные.

Запуск:
  python field_line_viz.py                       # демо-набор
  python field_line_viz.py SP_I1_1A 150,-8,-8 150,8,8   # своя линия P1→P2 (мм)
"""
import struct, os, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from field_compare_viz import fit_coeffs, ideal_field, peak_slice, ELEMS

HDR = 96
FIELD_DIR = os.path.join("workspace", "fields_bin")
FIGS = os.path.join("workspace", "bip_compare", "figs")
ORDERS = {dev: orders for dev, _t, orders, _n in ELEMS}
RU = {"QP_L_1A": "Квадруполь длинный", "QP_sh_1A": "Квадруполь короткий",
      "SP_I1_1A": "Секступоль (гл.)", "SP_I2_1A": "Секступоль (корр.)", "OP_1A": "Октуполь",
      "DQ_I1_1A": "Диполь-квадруполь I1", "DQ_I2_1A": "Диполь-квадруполь I2",
      "DQ_I3_1A": "Диполь-квадруполь I3", "DP_L_PM": "Диполь длинный", "DP_sh_PM": "Диполь короткий",
      "corrector_I1_1A": "Корректор I1", "corrector_I2_1A": "Корректор I2"}


def load(device):
    with open(os.path.join(FIELD_DIR, device + ".bin"), "rb") as f:
        b = f.read(HDR)
        mnx, mny, mnz, mxx, mxy, mxz, sx, sy, sz = struct.unpack_from("<9d", b, 0)
        nx, ny, nz = struct.unpack_from("<3i", b, 72)
        B = np.frombuffer(f.read(nx * ny * nz * 24), dtype="<f8").reshape(nx, ny, nz, 3)
    return dict(mnx=mnx, mny=mny, mnz=mnz, mxx=mxx, sx=sx, sy=sy, sz=sz, nx=nx, ny=ny, nz=nz), B


def trilin_sym(B, h, x, y, z):
    """Выборка с учётом YZ-симметрии карты: для x<mnx (за плоскостью x=0) поле зеркалится —
    Bx меняет знак, By/Bz чётны (как в FieldMap при use_symmetry=1)."""
    if x < h["mnx"]:
        bx, by, bz = trilin(B, h, 2 * h["mnx"] - x, y, z)   # отражение относительно x=mnx
        return np.array([-bx, by, bz])
    return trilin(B, h, x, y, z)


def trilin(B, h, x, y, z):
    """3D-трилинейная интерполяция вектора B в точке (x,y,z), м (только x≥mnx)."""
    fx = (x - h["mnx"]) / h["sx"]; fy = (y - h["mny"]) / h["sy"]; fz = (z - h["mnz"]) / h["sz"]
    i = int(np.clip(math.floor(fx), 0, h["nx"] - 2)); tx = fx - i
    j = int(np.clip(math.floor(fy), 0, h["ny"] - 2)); ty = fy - j
    k = int(np.clip(math.floor(fz), 0, h["nz"] - 2)); tz = fz - k
    out = np.zeros(3)
    for di, wx in ((0, 1 - tx), (1, tx)):
        for dj, wy in ((0, 1 - ty), (1, ty)):
            for dk, wz in ((0, 1 - tz), (1, tz)):
                out += B[i + di, j + dj, k + dk] * (wx * wy * wz)
    return out


def body_slice(h, B):
    """Срез ТЕЛА магнита — максимум ВНЕОСЕВого поперечного поля (а не осевого fringe-Bx)."""
    j0, k0 = h["ny"] // 2, h["nz"] // 2
    R = max(1, round(3e-3 / h["sy"]))
    reg = B[:, j0 - R:j0 + R + 1, k0 - R:k0 + R + 1, 1:3]      # By,Bz в кольце ±3мм
    mag = np.sqrt(np.mean(reg ** 2, axis=(1, 2, 3)))
    return int(np.argmax(mag)) if mag.max() > 0 else peak_slice(h, B)


def ideal_coeffs(device, h, B):
    """Коэффициенты идеального мультиполя, подогнанные в срезе ТЕЛА магнита (r≤5мм)."""
    orders = ORDERS.get(device, (2,))
    ix = body_slice(h, B)
    j0, k0 = h["ny"] // 2, h["nz"] // 2; R = 6
    js = slice(j0 - R, j0 + R + 1); ks = slice(k0 - R, k0 + R + 1)
    yg = (np.arange(h["ny"])[js] - j0) * h["sy"]; zg = (np.arange(h["nz"])[ks] - k0) * h["sz"]
    Zg, Yg = np.meshgrid(zg, yg)
    rr = np.sqrt(Yg ** 2 + Zg ** 2)
    return fit_coeffs(Yg, Zg, B[ix, js, ks, 1], B[ix, js, ks, 2], orders, rr <= 5e-3), orders


def field_line(device, P1, P2, n=240, path=None):
    """P1, P2 — (x,y,z) в мм (локальные координаты карты). Строит 3 панели Bx,By,Bz."""
    h, B = load(device)
    coeffs, orders = ideal_coeffs(device, h, B)
    j0, k0 = h["ny"] // 2, h["nz"] // 2
    yc = (j0) * h["sy"] + h["mny"]; zc = (k0) * h["sz"] + h["mnz"]   # координаты центра (м)
    p1 = np.array(P1) * 1e-3; p2 = np.array(P2) * 1e-3
    ts = np.linspace(0, 1, n)
    pts = np.array([p1 + t * (p2 - p1) for t in ts])
    s = np.linalg.norm(pts - p1, axis=1) * 1e3                       # длина вдоль линии, мм
    real = np.array([trilin(B, h, *p) for p in pts])                # реальное поле
    # идеальное: мультиполь по поперечным (y,z) относительно центра; Bx_ideal = 0
    ideal = np.zeros_like(real)
    for m, p in enumerate(pts):
        Y = np.array([[p[1] - yc]]); Z = np.array([[p[2] - zc]])
        iby, ibz = ideal_field(Y, Z, coeffs)
        ideal[m, 1] = iby[0, 0]; ideal[m, 2] = ibz[0, 0]
    fig, axes = plt.subplots(3, 1, figsize=(9.2, 8.4), sharex=True)
    names = ["Bx (продольная)", "By", "Bz"]
    for ci, (ax, nm) in enumerate(zip(axes, names)):
        ax.plot(s, real[:, ci] * 1e6, "-", color="tab:blue", lw=1.7, label="реальное (ANSYS)")
        ax.plot(s, ideal[:, ci] * 1e6, "--", color="tab:green", lw=1.6, label="идеальный мультиполь")
        ax.set_ylabel(nm + ", мкТл"); ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="best")
    axes[-1].set_xlabel("расстояние вдоль линии, мм")
    omod = "+".join(map(str, orders))
    fig.suptitle(f"{RU.get(device, device)} ({device}): срез поля вдоль линии\n"
                 f"P1=({P1[0]:g},{P1[1]:g},{P1[2]:g}) → P2=({P2[0]:g},{P2[1]:g},{P2[2]:g}) мм  "
                 f"(идеал = мультиполь N={omod}; Bx идеала ≡ 0)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    if path is None:
        path = os.path.join(FIGS, f"fig_line_{device}.png")
    fig.savefig(path, dpi=115); plt.close(fig)
    return path


# (подпись, смещение Y мм, смещение Z мм, тег файла)
PROFILE_OFFS = [("на оси (0,0)", 0.0, 0.0, "axis"),
                ("Y +2 мм", 2.0, 0.0, "Yp2"),
                ("Z +2 мм", 0.0, 2.0, "Zp2"),
                ("Y +8 мм", 8.0, 0.0, "Yp8"),
                ("Z +8 мм", 0.0, 8.0, "Zp8")]


def field_profile(device, y0mm, z0mm, lab, path):
    """Продольный профиль ОДНОГО смещения: три компоненты Bx,By,Bz вдоль X через плоскость
    симметрии (x от −maxX до +maxX, для x<0 — зеркало). Реал (сплошная) vs идеал (пунктир)."""
    h, B = load(device)
    coeffs, orders = ideal_coeffs(device, h, B)
    yc = h["mny"] + (h["ny"] // 2) * h["sy"]; zc = h["mnz"] + (h["nz"] // 2) * h["sz"]
    y0 = yc + y0mm * 1e-3; z0 = zc + z0mm * 1e-3
    xs = np.linspace(-h["mxx"], h["mxx"], min(2 * h["nx"], 800))   # СИММЕТРИЧНЫЙ диапазон через x=0
    real = np.array([trilin_sym(B, h, x, y0, z0) for x in xs])
    Y = np.array([[y0 - yc]]); Z = np.array([[z0 - zc]])
    iby, ibz = ideal_field(Y, Z, coeffs)
    ideal = [0.0, iby[0, 0], ibz[0, 0]]                            # Bx идеала ≡ 0
    fig, axes = plt.subplots(3, 1, figsize=(9.2, 8.0), sharex=True)
    comps = [("Bx (продольная)", "ОДД: −Bx(−x), на оси сим. = 0"), ("By", "чётна по x"),
             ("Bz", "чётна по x")]
    for ci, (ax, (nm, note)) in enumerate(zip(axes, comps)):
        ax.plot(xs * 1e3, real[:, ci] * 1e6, "-", color="tab:blue", lw=1.5, label="реальное (ANSYS)")
        ax.axhline(ideal[ci] * 1e6, ls="--", color="tab:green", lw=1.4, label="идеальный мультиполь")
        ax.axvline(0, color="0.5", ls=":", lw=1.2)                 # плоскость симметрии y-z
        ax.set_ylabel(nm + ", мкТл"); ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="best", title=note, title_fontsize=7)
    axes[-1].set_xlabel("X — продольная координата, мм  (пунктир x=0 — плоскость симметрии y-z)")
    omod = "+".join(map(str, orders))
    fig.suptitle(f"{RU.get(device, device)} ({device}): срез поля вдоль X, смещение {lab}\n"
                 f"диапазон через плоскость симметрии; идеал = мультиполь N={omod} (Bx ≡ 0)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=110); plt.close(fig)
    return path


ALL_DEVICES = ["QP_L_1A", "QP_sh_1A", "SP_I1_1A", "SP_I2_1A", "OP_1A",
               "DQ_I1_1A", "DQ_I2_1A", "DQ_I3_1A", "DP_L_PM", "DP_sh_PM",
               "corrector_I1_1A", "corrector_I2_1A"]


def main():
    os.makedirs(FIGS, exist_ok=True)
    if len(sys.argv) >= 2 and sys.argv[1] == "--profiles":
        for dev in ALL_DEVICES:
            for lab, y0, z0, tag in PROFILE_OFFS:
                p = field_profile(dev, y0, z0, lab,
                                  os.path.join(FIGS, f"fig_profile_{dev}_{tag}.png"))
                print("сохранено:", os.path.basename(p))
        return
    if len(sys.argv) >= 4:
        dev = sys.argv[1]
        P1 = [float(v) for v in sys.argv[2].split(",")]
        P2 = [float(v) for v in sys.argv[3].split(",")]
        print("сохранено:", field_line(dev, P1, P2))
        return
    def body_x_mm(dev):
        h, B = load(dev)
        return (body_slice(h, B) * h["sx"] + h["mnx"]) * 1e3      # x тела магнита, мм
    # демо: поперечная диагональ в срезе ТЕЛА + продольная линия (fringe)
    xs = body_x_mm("SP_I1_1A"); xq = body_x_mm("QP_L_1A"); xd = body_x_mm("DQ_I1_1A")
    demos = [
        ("SP_I1_1A", [xs, -8, -8], [xs, 8, 8], "_diag"),       # секступоль: диагональ в теле
        ("QP_L_1A", [xq, -8, -8], [xq, 8, 8], "_diag"),        # квадруполь: диагональ в теле
        ("QP_L_1A", [0, 5, 0], [400, 5, 0], "_long"),          # квадруполь: продольная (fringe, Bx)
        ("DQ_I1_1A", [xd, -8, 0], [xd, 8, 0], "_diag"),        # DQ: поперечная по Y в теле
    ]
    for dev, P1, P2, tag in demos:
        p = field_line(dev, P1, P2, path=os.path.join(FIGS, f"fig_line_{dev}{tag}.png"))
        print(f"сохранено: {os.path.basename(p)}  (тело x={P1[0]:.0f}мм)" if tag == "_diag"
              else f"сохранено: {os.path.basename(p)}")


if __name__ == "__main__":
    main()
