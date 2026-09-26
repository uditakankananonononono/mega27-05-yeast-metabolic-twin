#!/usr/bin/env python3
"""Amendment-2 sensitivity sweep: one-at-a-time ASSUMED parameters.

7 parameter sets (nominal + 2 per swept parameter), each on its own A1
grid curves and the locked 1,000-environment subsample (every 5th of the
5,000 LHS). Checkpointed: rows append to sweep_<name>.csv, reruns skip
finished envs. --budget seconds bounds each invocation.
"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from yeasttwin.model import load_model  # noqa: E402
from yeasttwin.climate.apply import prepare  # noqa: E402
from yeasttwin.climate.environment import Environment, grid_axes, lhs_sample  # noqa: E402
from yeasttwin.climate.objectives import ethanol_yield  # noqa: E402
from yeasttwin.climate.parameters import (EthanolModel, OsmoticModel,  # noqa: E402
                                          StressParameters, TemperatureModel)

SETS = {
    "nominal": StressParameters(),
    "kgly_lo": StressParameters(osmotic=OsmoticModel(k_gly=0.5)),
    "kgly_hi": StressParameters(osmotic=OsmoticModel(k_gly=2.0)),
    "n_lo": StressParameters(ethanol=EthanolModel(n=0.75)),
    "n_hi": StressParameters(ethanol=EthanolModel(n=2.25)),
    "tmax_lo": StressParameters(temperature=TemperatureModel(t_max=40.0)),
    "tmax_hi": StressParameters(temperature=TemperatureModel(t_max=44.0)),
}

AXES = ("temperature_c", "ethanol_pct", "osmotic_m", "nitrogen_frac")


def a1_curves(model, ref, params):
    curves = {}
    for axis in AXES:
        pts = []
        for v in grid_axes()[axis]:
            env = Environment(**{axis: float(v)})
            r = ethanol_yield(model, env, ref, params)
            pts.append((float(v), r.ethanol))
        curves[axis] = pts
    ref_e = curves["temperature_c"][2][1]  # 30 C is index 2
    if ref_e <= 0:
        raise RuntimeError("zero reference in sweep")
    return {a: [(x, y / ref_e) for x, y in pts] for a, pts in curves.items()}, ref_e


def interp(pts, x):
    if x <= pts[0][0]:
        return pts[0][1]
    if x >= pts[-1][0]:
        return pts[-1][1]
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= x <= x1:
            return y0 + (x - x0) * (y1 - y0) / (x1 - x0)
    return pts[-1][1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=sorted(SETS))
    ap.add_argument("--budget", type=float, default=95.0)
    args = ap.parse_args()
    params = SETS[args.set]
    t0 = time.time()

    envs = lhs_sample()[::5]  # locked 1,000-env subsample
    out = ROOT / "results" / "climate" / f"sweep_{args.set}.csv"
    done = set()
    if out.exists():
        for r in csv.DictReader(open(out)):
            done.add(int(r["idx"]))

    model = load_model()
    ref = prepare(model)
    curves, ref_e = a1_curves(model, ref, params)

    mode = "a" if out.exists() else "w"
    with open(out, mode, newline="") as fh:
        w = csv.writer(fh)
        if mode == "w":
            w.writerow(["idx", "relative_yield", "collapsed", "deviation"])
        for i, env in enumerate(envs):
            if i in done:
                continue
            if time.time() - t0 > args.budget:
                print(f"BUDGET-EXIT set={args.set} done={len(done)}")
                return
            r = ethanol_yield(model, env, ref, params)
            rel = r.ethanol / ref_e if r.feasible else 0.0
            indep = 1.0
            for a in AXES:
                indep *= interp(curves[a], getattr(env, a))
            w.writerow([i, rel, (not r.feasible) or rel < 0.2, rel - indep])
            done.add(i)
    # complete: summarize
    rows = list(csv.DictReader(open(out)))
    collapsed = [r for r in rows if r["collapsed"] == "True"]
    neg = [r for r in collapsed if float(r["deviation"]) < 0]
    summary = {
        "set": args.set,
        "n_envs": len(rows),
        "n_collapsed": len(collapsed),
        "collapsed_fraction": len(collapsed) / len(rows),
        "fraction_negative_deviation_in_collapsed": (len(neg) / len(collapsed)) if collapsed else None,
        "mean_deviation_collapsed": (sum(float(r["deviation"]) for r in collapsed) / len(collapsed)) if collapsed else None,
    }
    (ROOT / "results" / "climate" / f"sweep_{args.set}_summary.json").write_text(json.dumps(summary, indent=2))
    print("DONE", json.dumps(summary))


if __name__ == "__main__":
    main()
