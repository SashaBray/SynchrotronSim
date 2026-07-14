# -*- coding: utf-8 -*-
"""Генерация СИНТЕТИЧЕСКИХ активов для новых элементов библиотеки (флаг --gen-assets):
  gen_field_map — карта поля масштабированием карты эталона (ГИБРИДНЫЙ ресайз: тело ×c_body под нужную
    L_eff, фриндж сохранён — из исследования масштабирования карт). Выход — компактный .bin.
  gen_stl       — 3D-модель продольным масштабом STL эталона (растяжение/сжатие вдоль длинной оси).
Обе функции устойчивы: при любой ошибке возвращают None (актив не создаётся, поле остаётся '-')."""
import os
import struct
import glob
import numpy as np
from field_quality import load
import _qp_scale as qs
import _qp_scale_hybrid as hy

GEO = os.path.join("workspace", "geometry")
FB = os.path.join("workspace", "fields_bin")
FR = os.path.join("workspace", "fields_raw")


HDR = 96                                                       # размер шапки .bin (как field_quality.HDR)


def _write_bin(path, h, B, nominal_arg=1.0):
    """Запись карты в .bin формата C++ GridHeader (AnsysConverter.h), 96 байт:
      <9d> minX,minY,minZ,maxX,maxY,maxZ,stepX,stepY,stepZ | <3i> nx,ny,nz | <4x> паддинг | <d> nominalArg.
    КРИТИЧНО: nominalArg (ток снятия карты) на смещении 88 — C++ делит на него (currentScale=arg_val/
    nominalArg); при 0 получается NaN. Для 1А-карт nominal_arg=1.0 (как у базовых)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    mxy = h["mny"] + (h["ny"] - 1) * h["sy"]; mxz = h["mnz"] + (h["nz"] - 1) * h["sz"]
    hdr = (struct.pack("<9d", h["mnx"], h["mny"], h["mnz"], h["mxx"], mxy, mxz, h["sx"], h["sy"], h["sz"])
           + struct.pack("<3i", int(h["nx"]), int(h["ny"]), int(h["nz"])))       # 84 байта
    hdr += b"\x00" * (88 - len(hdr))                           # выравнивание double: паддинг до смещения 88
    hdr += struct.pack("<d", float(nominal_arg))              # nominalArg на смещении 88 (как в C++ struct)
    with open(path, "wb") as f:
        f.write(hdr); f.write(np.ascontiguousarray(B, dtype="<f8").tobytes())


def gen_field_map(ref, nid, new_eff, ref_row, col):
    """Карта поля синтетического элемента: гибридный ресайз тела карты эталона под L_eff=new_eff.
    Возврат: путь .bin (относит. fields_bin) для map_bin_path, либо None."""
    try:
        raw = ref_row[col["map_raw_path"]]
        src = os.path.join(FR, os.path.basename(raw)) if raw not in ("-", "") else None
        if not src or not os.path.exists(src):
            return None
        h, B = load(src)
        if h["mnx"] >= -1e-9:                                  # полукарта X≥0 → полная (зеркало)
            h, B = qs.mirror_full(h, B)
        xg = h["mnx"] + np.arange(h["nx"]) * h["sx"]
        yg = h["mny"] + np.arange(h["ny"]) * h["sy"]; zg = h["mnz"] + np.arange(h["nz"]) * h["sz"]
        kc = np.argmin(np.abs(zg)); jp = np.argmin(np.abs(yg - 0.35 * max(abs(yg[0]), yg[-1])))
        prof = np.sqrt((B[:, jp, kc, :] ** 2).sum(-1))         # |B| в приосевой точке: плоско в теле, спад во фриндже
        Gm = prof.max()
        if Gm <= 0:
            return None
        pos = np.where(xg >= 0)[0]; xp = xg[pos]; pp = prof[pos] / Gm
        Xb = next((xp[i] for i in range(len(pp)) if pp[i] < 0.9), xp[-1])    # полуширина тела по профилю
        phys_ref = float(ref_row[col["phys_length"]])
        if Xb < 0.4 * phys_ref:                                 # профиль не выделил тело (combined-func) → по металлу
            Xb = 0.5 * phys_ref
        L_body = 2 * Xb
        # ДОБАВИТЬ/УБРАТЬ ровно ΔL=new_eff−eff_ref тела: истинная L_eff меняется на ΔL (линейно),
        # независимо от меры профиля. c_fringe=1 (фриндж сохранён).
        eff_ref = float(ref_row[col["eff_length"]]); dL = new_eff - eff_ref
        cb = 1.0 + dL / L_body if L_body > 1e-6 else 1.0
        cb = min(max(cb, 0.05), 20.0)
        h2, B2 = hy.resize(h, B, cb, 1.0, Xb)                  # тело ×c_body, фриндж сохранён
        # ХРАНИМ ПОЛУКАРТОЙ X≥0 (как эталоны: use_symmetry=1). Карта симметрична относительно X=0
        # (зеркалили полукарту → ресайз тела симметричен), поэтому C++ восстановит −половину зеркалом.
        # Вдвое меньше диска и точное соответствие соглашению библиотеки.
        xg2 = h2["mnx"] + np.arange(h2["nx"]) * h2["sx"]
        k0 = int(np.searchsorted(xg2, -1e-9))                  # первый индекс X≥0
        if 0 < k0 < h2["nx"]:
            B2 = np.ascontiguousarray(B2[k0:, :, :, :])
            h2 = dict(h2); h2["nx"] = B2.shape[0]
            h2["mnx"] = float(xg2[k0]); h2["mxx"] = float(xg2[-1])
        na = float(ref_row[col["nominal_arg"]]) if "nominal_arg" in col and ref_row[col["nominal_arg"]] \
            not in ("", "-") else 1.0                          # ток снятия карты эталона (для currentScale)
        rel = f"synthetic/{nid}.bin"
        _write_bin(os.path.join(FB, "synthetic", f"{nid}.bin"), h2, B2, na)
        return rel
    except Exception as e:
        print(f"      [gen_field_map {nid}] пропуск: {type(e).__name__}: {e}")
        return None


def gen_stl(ref, nid, new_eff, ref_row, col):
    """3D-модель синтетического элемента: STL эталона, растянутый/сжатый вдоль длинной оси под L_phys.
    Возврат: путь .STL (относит. geometry) для stl_path, либо None."""
    try:
        stl = ref_row[col["stl_path"]]
        if stl in ("-", ""):
            return None
        src = os.path.join(GEO, stl)
        if not os.path.exists(src):                            # регистр (.stl в БД vs .STL на диске)
            cand = glob.glob(os.path.join(GEO, os.path.splitext(stl)[0] + ".*"))
            src = cand[0] if cand else None
        if not src or not os.path.exists(src):
            return None
        ref_eff = float(ref_row[col["eff_length"]]); ref_phys = float(ref_row[col["phys_length"]])
        new_phys = ref_phys + (new_eff - ref_eff)              # сплайс тела (как в библиотеке)
        ratio = new_phys / ref_phys if ref_phys > 1e-6 else 1.0
        with open(src, "rb") as f:
            head = f.read(80); n = struct.unpack("<I", f.read(4))[0]; raw = f.read(n * 50)
        dt = np.dtype([("n", "<f4", 3), ("v1", "<f4", 3), ("v2", "<f4", 3), ("v3", "<f4", 3), ("a", "<u2")])
        tris = np.frombuffer(raw, dtype=dt, count=n).copy()
        V = np.concatenate([tris["v1"], tris["v2"], tris["v3"]])
        ax = 0                                                 # ПРОДОЛЬНАЯ (пучковая) ось STL = X, как у полевых
        # карт (из той же КД): у всех магнитов X-размах ≈ длине, Y/Z — поперечные (ярмо+катушки). НЕ argmax bbox.
        ctr = 0.5 * (V.max(0)[ax] + V.min(0)[ax])
        for key in ("v1", "v2", "v3"):
            arr = tris[key].copy()                             # поле как (n,3); правим и присваиваем ЦЕЛИКОМ
            arr[:, ax] = ctr + (arr[:, ax] - ctr) * ratio      # масштаб вдоль длинной оси от центра
            tris[key] = arr
        e1 = tris["v2"] - tris["v1"]; e2 = tris["v3"] - tris["v1"]; nrm = np.cross(e1, e2)
        ln = np.linalg.norm(nrm, axis=1, keepdims=True)
        tris["n"] = np.divide(nrm, ln, out=np.zeros_like(nrm), where=ln > 0)   # безопасно (вырожд. треуг.)
        out = os.path.join(GEO, "synthetic", f"{nid}.STL"); os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "wb") as f:
            f.write(head); f.write(struct.pack("<I", n)); f.write(tris.tobytes())
        return f"synthetic/{nid}.STL"
    except Exception as e:
        print(f"      [gen_stl {nid}] пропуск: {type(e).__name__}: {e}")
        return None
