# -*- coding: utf-8 -*-
"""Очистка синтетических карт ДИПОЛЕЙ от паразитного медианного поперечного поля.

У медианно-симметричного магнита (медиана = плоскость кольца, z=0 в карте) продольное Bx и
поперечно-горизонтальное By должны быть НЕЧЁТНЫ по z (=0 на медиане), а вертикальное Bz — ЧЁТНО.
Синтетические карты диполей DP несут нефизичный DC-сдвиг By≈+0.11мТл, Bx≈−0.18мТл на медиане
(артефакт эталонного ANSYS-решения, скопированный ресайзом) — он даёт систематический ВЕРТИКАЛЬНЫЙ
кик каждому обороту (By/Bz≈3.6e−4 × полный поворот 2π ⇒ ~2.3мрад увода орбиты, пучок гибнет ~1.4 SP).

Фикс: наложить точную медианную симметрию — Bx,By ← нечётная часть по z, Bz ← чётная часть.
Заголовок .bin (96 байт, вкл. nominalArg на смещении 88) сохраняется бит-в-бит. Бэкап оригинала.

Запуск: python clean_dipole_maps.py            # чистит все DP_*.bin в synthetic
        python clean_dipole_maps.py <dev> ...  # конкретные device_type_id
"""
import os
import sys
import glob
import shutil
import struct
import numpy as np

FB = os.path.join("workspace", "fields_bin", "synthetic")
HDR = 96


def _read(path):
    with open(path, "rb") as f:
        head = f.read(HDR)
        mnx, mny, mnz, mxx, mxy, mxz, sx, sy, sz = struct.unpack_from("<9d", head, 0)
        nx, ny, nz = struct.unpack_from("<3i", head, 72)
        B = np.frombuffer(f.read(nx * ny * nz * 24), dtype="<f8").reshape(nx, ny, nz, 3).copy()
    return head, dict(mnz=mnz, mxz=mxz, sz=sz, nx=nx, ny=ny, nz=nz), B


def clean(path):
    head, h, B = _read(path)
    nz = h["nz"]
    if abs(h["mnz"] + h["mxz"]) > 1e-6:                       # z-сетка обязана быть симметрична (−mxz..mxz)
        print(f"  ПРОПУСК {os.path.basename(path)}: z-сетка не симметрична ({h['mnz']:.4f}..{h['mxz']:.4f})")
        return False
    ix, jy, kz = h["nx"] // 2, h["ny"] // 2, nz // 2
    by0, bx0 = B[ix, jy, kz, 1] * 1e3, B[ix, jy, kz, 0] * 1e3
    Bflip = B[:, :, ::-1, :]                                  # отражение по z (сетка симм. ⇒ индексы совпадают)
    Bc = B.copy()
    Bc[..., 0] = 0.5 * (B[..., 0] - Bflip[..., 0])            # Bx ← нечётная часть
    Bc[..., 1] = 0.5 * (B[..., 1] - Bflip[..., 1])            # By ← нечётная часть
    Bc[..., 2] = 0.5 * (B[..., 2] + Bflip[..., 2])            # Bz ← чётная часть
    by1, bx1 = Bc[ix, jy, kz, 1] * 1e3, Bc[ix, jy, kz, 0] * 1e3
    bak = path + ".orig"
    if not os.path.exists(bak):
        shutil.copy(path, bak)
    with open(path, "wb") as f:
        f.write(head)
        f.write(np.ascontiguousarray(Bc, dtype="<f8").tobytes())
    print(f"  {os.path.basename(path):26} медиана By {by0:+.4f}→{by1:+.4f}  Bx {bx0:+.4f}→{bx1:+.4f} мТл  [бэкап .orig]")
    return True


def main(argv):
    if argv:
        paths = [os.path.join(FB, a if a.endswith(".bin") else a + ".bin") for a in argv]
    else:
        paths = sorted(glob.glob(os.path.join(FB, "DP_*.bin")))
    print(f"Очистка {len(paths)} карт диполей (медианная симметрия Bx,By нечёт / Bz чёт по z):")
    n = sum(clean(p) for p in paths if os.path.exists(p))
    print(f"[OK] очищено {n} карт.")


if __name__ == "__main__":
    main(sys.argv[1:])
