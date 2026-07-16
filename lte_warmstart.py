# -*- coding: utf-8 -*-
"""Warm-start для матричного оптимизатора: проектные силы из .lte → физкоэффициенты магнитов.

Идея: .lte SYLA_desc задаёт K1/K2/K3/ANGLE каждого именованного магнита, а SYLA_line — порядок
магнитов в стандартном супериоде. Наши сгруппированные магниты (lattice_configs 105, со-локованные
обмотки слиты) идут в ТОМ ЖЕ порядке. Сопоставляем проходом с проверкой типа (kind-matched walk):
где типы совпадают — берём проектную силу; спец-случаи (SH-секступоль .lte ↔ наш корректор;
junction-комбайнер DMB/SBM) помечаем и оставляем на fallback (оптимизатор дотюнит).

Конвенция сил (elegant KQUAD/KSEXT/KOCT == наши k):  k1=K1, k2=K2, k3=K3  ⇒  физкоэф.:
   G = K1/QP,  S = K2/QP,  O = K3/QP    (QP = q/p₀ = bm.QP).  θ_bend = ANGLE (геометрия, у диполей/DQ).
"""
import os, re, csv
import bip_matrix as bm

QP = bm.QP
LTE = os.path.join("workspace", "lte")
TB = os.path.join("workspace", "tables")
_DESC_NAME = "SYLA_desc.lte"
_LINE_NAME = "SYLA_line.lte"


def set_lte_dir(path):
    """Переключить источник .lte на папку конфигурации (ищет *_desc.lte/*_line.lte). Без копий в корень."""
    global LTE, _DESC_NAME, _LINE_NAME
    import glob as _glob
    LTE = path
    for pat, var in (("*_desc.lte", "_DESC_NAME"), ("*_line.lte", "_LINE_NAME")):
        hits = sorted(_glob.glob(os.path.join(path, pat)))
        if hits:
            globals()[var] = os.path.basename(hits[0])


# ---------------- разбор .lte ----------------
def parse_desc():
    """SYLA_desc.lte → {name: {kind, K1, K2, K3, angle, L}}. kind ∈ quad|sext|oct|dq|dipole."""
    txt = open(os.path.join(LTE, _DESC_NAME), encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"^([A-Za-z0-9_]+)\s*:\s*([A-Za-z_]+)\s*,(.*)$", txt, re.M):
        name, typ, rest = m.group(1), m.group(2).upper(), m.group(3)
        kv = {}
        for pair in rest.split(","):
            if "=" in pair:
                k, v = pair.split("=", 1)
                try: kv[k.strip().upper()] = float(v)
                except ValueError: pass
        K1, K2, K3 = kv.get("K1", 0.0), kv.get("K2", 0.0), kv.get("K3", 0.0)
        ang, L = kv.get("ANGLE", 0.0), kv.get("L", 0.0)
        if typ == "KQUAD": kind = "quad"
        elif typ == "KSEXT": kind = "sext"
        elif typ == "KOCT": kind = "oct"
        elif typ in ("CSBEND", "SBEND", "RBEN"):
            kind = "dq" if abs(K1) > 1e-12 else "dipole"
        else: kind = "other"
        out[name] = dict(kind=kind, K1=K1, K2=K2, K3=K3, angle=ang, L=L, typ=typ)
    return out


def parse_line_superperiod():
    """SYLA_line.lte → список имён магнитов ОДНОГО стандартного супериода (до первого повтора QF1A)."""
    txt = open(os.path.join(LTE, _LINE_NAME), encoding="utf-8").read()
    names = [t.strip() for t in txt.replace("\n", ",").split(",") if t.strip()]
    # супериод = от первого QF1A до следующего QF1A (не включая)
    i0 = names.index("QF1A")
    i1 = names.index("QF1A", i0 + 1)
    seq = names[i0:i1]
    mags = [n for n in seq if not (n.startswith("DR_") or n == "BPM" or n.startswith("MARK"))]
    return mags


def parse_line_full():
    """SYLA_line.lte → список ВСЕХ магнитов ПОЛНОГО кольца по порядку (38 SP + junction)."""
    txt = open(os.path.join(LTE, _LINE_NAME), encoding="utf-8").read()
    names = [t.strip() for t in txt.replace("\n", ",").split(",") if t.strip()]
    return [n for n in names if not (n.startswith("DR_") or n == "BPM" or n.startswith("MARK"))]


# ---------------- наши сгруппированные магниты ----------------
def _load_dev_kind():
    """kind по устройству из element_transforms.csv (включая СИНТЕТИЧЕСКИЕ элементы); kicker→corrector."""
    kmap = {"kicker": "corrector"}
    out = {}
    for r in csv.DictReader(open(os.path.join(TB, "element_transforms.csv"), encoding="utf-8"), delimiter=";"):
        out[r["device_type_id"]] = kmap.get(r["kind"], r["kind"])
    return out


DEV_KIND = _load_dev_kind()


def grouped_magnets(lattice_id=105):
    """Инстансы решётки → список групп со-локованных обмоток (как build_line трекера)."""
    inst = [r for r in csv.DictReader(open(os.path.join(TB, "lattice_configs.csv"), encoding="utf-8"),
                                      delimiter=";") if r["lattice_id"] == str(lattice_id)]
    inst.sort(key=lambda r: int(r["instance_id"]))
    grp = []; i = 0
    while i < len(inst):
        r = inst[i]; g = [r]
        while i + 1 < len(inst) and (inst[i + 1]["x"], inst[i + 1]["y"], inst[i + 1]["z"]) == (r["x"], r["y"], r["z"]):
            g.append(inst[i + 1]); i += 1
        grp.append(g); i += 1
    return grp


def our_kind(group):
    return DEV_KIND.get(group[0]["device_type_id"], "other")


def load_eff():
    """Эффективные длины (полевой интеграл) по типу устройства — для сохранения проектного интеграла."""
    return {r["device_type_id"]: float(r["eff_length"])
            for r in csv.DictReader(open(os.path.join(TB, "hardware_library.csv"), encoding="utf-8"),
                                    delimiter=";")}


# ---------------- сопоставление (kind-matched walk) ----------------
def _walk(lte_seq, grp, desc, eff, window=5):
    """Проход: сопоставляет каждой нашей группе магнит .lte совместимого типа (в порядке), берёт
    проектный K с сохранением интеграла K·Leff=K1·L_lte. corrector ↔ SH-секступоль (K=0). Возвращает
    список dict(idx, devs, kind, lte, G, S, O, angle, note)."""
    lte_kinds = [desc.get(n, {}).get("kind", "other") for n in lte_seq]
    res = []
    j = 0
    for i, g in enumerate(grp):
        ok = our_kind(g)
        devs = "+".join(x["device_type_id"] for x in g)
        note = ""
        if ok == "corrector":
            if j < len(lte_kinds) and lte_kinds[j] == "sext":
                note = f"корректор ↔ .lte секступоль {lte_seq[j]}"; j += 1
            else:
                note = "корректор (нет .lte)"
            res.append(dict(idx=i, devs=devs, kind=ok, lte=None, G=0.0, S=0.0, O=0.0, angle=0.0, note=note))
            continue
        matched = -1
        for jj in range(j, min(j + window, len(lte_kinds))):
            lk = lte_kinds[jj]
            if lk == ok or (ok == "dipole" and lk in ("dipole", "dq")) or (ok == "dq" and lk in ("dq", "dipole")):
                matched = jj; break
        if matched >= 0:
            name = lte_seq[matched]; d = desc[name]
            if matched > j:
                note = f"пропущено: {lte_seq[j:matched]}"
            # СОХРАНЕНИЕ ИНТЕГРАЛА: k·Leff = K·L_lte ⇒ coeff = (K/QP)·(L_lte/Leff).
            Leff = eff.get(g[0]["device_type_id"], d["L"] or 1.0)
            scale = (d["L"] / Leff) if (Leff and d["L"]) else 1.0
            G = (d["K1"] / QP) * scale if ok in ("quad", "dq") else 0.0
            S = (d["K2"] / QP) * scale if ok in ("sext", "dq") else 0.0
            O = (d["K3"] / QP) * scale if ok == "oct" else 0.0
            res.append(dict(idx=i, devs=devs, kind=ok, lte=name, G=G, S=S, O=O, angle=d["angle"],
                            note=note, Llte=d["L"], Leff=Leff, scale=scale))
            j = matched + 1
        else:
            note = f"НЕ НАЙДЕНО (j={j}, next={lte_seq[j] if j < len(lte_seq) else '—'})"
            res.append(dict(idx=i, devs=devs, kind=ok, lte=None, G=0.0, S=0.0, O=0.0, angle=0.0, note=note))
    return res


def build_mapping(lattice_id=105, verbose=True):
    """Warm-start одного стандартного супериода (реш. 105): проектные K → физкоэф. на 52 группы."""
    desc = parse_desc(); eff = load_eff()
    lte_seq = parse_line_superperiod()
    grp = grouped_magnets(lattice_id)
    res = _walk(lte_seq, grp, desc, eff)

    if verbose:
        print(f"QP=q/p₀={QP:.4f}. Групп: {len(grp)}, магнитов .lte: {len(lte_seq)}. Сопоставление:")
        print(f"{'i':>3} {'наши устройства':30} {'тип':9} {'.lte':10} {'G':>9} {'S':>9} {'O':>10} note")
        for r in res:
            print(f"{r['idx']:3} {r['devs']:30} {r['kind']:9} {str(r['lte']):10} "
                  f"{r['G']:+9.3f} {r['S']:+9.2f} {r['O']:+10.1f} {r['note']}")
        used = set(r["lte"] for r in res if r["lte"])
        miss = [n for n in lte_seq if n not in used and desc.get(n, {}).get("kind") in ("quad", "sext", "oct", "dq")]
        print(f"\nНе привязанные силовые магниты .lte: {miss}")
    return res


def build_mapping_full(lattice_id=101, verbose=True):
    """Warm-start ВСЕГО кольца (реш. 101, 38 SP + junction): проектные K из полной линии .lte на каждую
    сгруппированную группу магнитов. Возвращает список dict по всем группам кольца (в порядке)."""
    desc = parse_desc(); eff = load_eff()
    lte_seq = parse_line_full()
    grp = grouped_magnets(lattice_id)
    res = _walk(lte_seq, grp, desc, eff)
    if verbose:
        from collections import Counter
        nmatched = sum(1 for r in res if r["lte"])
        by_kind = Counter(r["kind"] for r in res)
        nznote = [r for r in res if r["note"].startswith("НЕ НАЙДЕНО")]
        print(f"QP={QP:.4f}. Кольцо {lattice_id}: групп {len(grp)}, магнитов .lte {len(lte_seq)}. "
              f"Привязано {nmatched}/{len(grp)}.")
        print(f"  по типам групп: {dict(by_kind)}")
        print(f"  НЕ найдено соответствий: {len(nznote)}")
        for r in nznote[:10]:
            print(f"    idx {r['idx']} {r['devs']}: {r['note']}")
    return res


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "full":
        build_mapping_full()
    else:
        build_mapping()
