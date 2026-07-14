#!/usr/bin/env python3
"""tune.py — одиночная оптимизация с РУБЕЖОМ, заданным аргументом (Nsp), БЕЗ плодения конфигов.

Финиш-плоскость ставится на gid = round(N*69) кольца (lattice 107), loss по магнитам
SP1..ceil(N). Та же машинерия, что у continuation_tuner (решётка/веса/группы), но без
авто-цикла — один проход с явным рубежом. Веса по умолчанию 1/1/10 (position/parallelism/
envelope_decay), как в текущих сериях.

Примеры:
  python tune.py --finish 1sp   --budget 5000
  python tune.py --finish 2.5sp --budget 5000
  python tune.py --finish 1sp   --budget 5000 --search monte_carlo --kick-every 500 --kick-scale 1.0
  python tune.py --finish 3sp   --budget 8000 --svd-angle-weight 3
"""
import argparse, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import continuation_tuner as ct


def parse_finish(s: str) -> float:
    m = re.fullmatch(r"\s*([0-9]*\.?[0-9]+)\s*sp\s*", s, re.I)
    if not m:
        raise argparse.ArgumentTypeError(f"--finish ожидает формат 'Nsp' (напр. 1sp, 2.5sp), получено '{s}'")
    return float(m.group(1))


def main() -> None:
    ap = argparse.ArgumentParser(description="Оптимизация с рубежом-аргументом (Nsp), без плодения конфигов.")
    ap.add_argument("--finish", type=parse_finish, required=True,
                    help="РУБЕЖ в супериодах: '1sp', '2.5sp'. Финиш-плоскость на gid=round(N*69).")
    ap.add_argument("--budget", type=int, required=True, help="бюджет симуляций")
    ap.add_argument("--search", choices=["line_search", "monte_carlo", "svd"], default="line_search",
                    help="метод глобального поиска")
    ap.add_argument("--parallelism-weight", type=float, default=1.0)
    ap.add_argument("--envelope-weight", type=float, default=10.0)
    ap.add_argument("--svd-angle-weight", type=float, default=0.0)
    ap.add_argument("--ls-max-rounds", type=int, default=12, help="line_search: потолок раундов дозума")
    ap.add_argument("--kick-every", type=int, default=0, help="monte_carlo: кик каждые N итераций")
    ap.add_argument("--kick-scale", type=float, default=1.0, help="monte_carlo: масштаб кика")
    ap.add_argument("--check-reach", action="store_true",
                    help="после оптимизации прогнать всё кольцо и показать достигнутый рубеж (медленно)")
    args = ap.parse_args()

    endpoints = ct.load_endpoints()
    raw_gid = max(1, min(ct.RING_END_GID, round(args.finish * ct.SP_PERIOD)))
    finish_gid = ct.snap_gid(endpoints, raw_gid)        # к концу co-located группы (не резать DQ пополам)
    n_sp = min(10, (finish_gid - 1) // ct.SP_PERIOD + 1)
    snap_note = f" (снап {raw_gid}->{finish_gid})" if finish_gid != raw_gid else ""
    print(f">>> РУБЕЖ {args.finish} SP = {ct._gid_lbl(finish_gid)}{snap_note}; loss по SP1..{n_sp}; "
          f"бюджет {args.budget}; метод {args.search}; веса pos=1/par={args.parallelism_weight}/"
          f"env={args.envelope_weight}" + (f"/svd_angle={args.svd_angle_weight}" if args.svd_angle_weight else ""),
          flush=True)

    ct.run_polish(endpoints, finish_gid, n_sp, args.budget, args.search,
                  kick_every=args.kick_every if args.search == "monte_carlo" else 0,
                  kick_scale=args.kick_scale,
                  ls_max_rounds=args.ls_max_rounds,
                  env_w=args.envelope_weight, par_w=args.parallelism_weight,
                  svd_angle_w=args.svd_angle_weight)

    if args.check_reach:
        reach, mx = ct._assemble_track_reach(ct.load_endpoints(), 0.012)
        print(f">>> ПОСЛЕ: пучок доходит до {ct._gid_lbl(reach)}, max откл. по кольцу {mx:.1f} мм", flush=True)


if __name__ == "__main__":
    main()
