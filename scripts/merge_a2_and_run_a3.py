#!/usr/bin/env python3
"""Merge A2 slices, then A3: Bliss-independence non-additivity (pre-reg 7.A3).

Independence prediction per environment = product of the four single-stress
relative yields (A1 curves, linearly interpolated to the LHS values).
Deviation = simulated relative yield - independence prediction. Locked
metric: fraction of collapsed-region environments with deviation < 0, with
a bootstrap 95% CI over environments. The pre-registered parameter-
uncertainty bootstrap (ASSUMED-parameter sweeps) is a separate later run.
"""
import csv
import glob
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "climate"

# --- merge ---
rows = []
for f in sorted(glob.glob(str(OUT / "a2_part_*.csv"))):
    rows.extend(csv.DictReader(open(f)))
seen = {}
for r in rows:
    seen[int(r["idx"])] = r
assert len(seen) == 5000, f"coverage gap: {len(seen)} != 5000"
merged = [seen[i] for i in range(5000)]
with open(OUT / "a2_collapse_surface.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(merged[0].keys()))
    w.writeheader()
    w.writerows(merged)

# --- A1 curves ---
a1 = list(csv.DictReader(open(OUT / "a1_dose_response.csv")))
curves = {}
for axis in ("temperature_c", "ethanol_pct", "osmotic_m", "nitrogen_frac"):
    pts = sorted(((float(r["value"]), float(r["relative_yield"]))
                  for r in a1 if r["axis"] == axis))
    curves[axis] = pts


def interp(axis, x):
    pts = curves[axis]
    if x <= pts[0][0]:
        return pts[0][1]
    if x >= pts[-1][0]:
        return pts[-1][1]
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= x <= x1:
            t = (x - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)
    return pts[-1][1]


# --- A3 deviations ---
out_rows = []
for r in merged:
    rel = float(r["relative_yield"])
    indep = (interp("temperature_c", float(r["temperature_c"]))
             * interp("ethanol_pct", float(r["ethanol_pct"]))
             * interp("osmotic_m", float(r["osmotic_m"]))
             * interp("nitrogen_frac", float(r["nitrogen_frac"])))
    collapsed = r["collapsed"] == "True"
    out_rows.append({"idx": r["idx"], "relative_yield": rel,
                     "independence_pred": indep, "deviation": rel - indep,
                     "collapsed": collapsed})
with open(OUT / "a3_nonadditivity.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    w.writerows(out_rows)

collapsed = [r for r in out_rows if r["collapsed"]]
neg = [r for r in collapsed if r["deviation"] < 0]
rng = random.Random(20260926)
n = len(collapsed)
boot = []
for _ in range(2000):
    samp = [collapsed[rng.randrange(n)] for _ in range(n)] if n else []
    boot.append(sum(1 for r in samp if r["deviation"] < 0) / n if n else 0.0)
boot.sort()
summary = {
    "analysis": "A3 Bliss-independence non-additivity (WT, nominal parameters)",
    "prereg_sha256_amended": "5256070ccb71617a6f88b15550f392332c2ea342e0e2c1e04782a0334a1bd8ae",
    "n_environments": len(out_rows),
    "n_collapsed": n,
    "collapsed_fraction": n / len(out_rows),
    "collapsed_with_negative_deviation": len(neg),
    "fraction_negative_deviation_in_collapsed": (len(neg) / n) if n else None,
    "bootstrap95_ci_fraction_negative": [boot[50], boot[1949]] if n else None,
    "mean_deviation_collapsed": (sum(r["deviation"] for r in collapsed) / n) if n else None,
    "parameter_uncertainty_bootstrap": "PENDING (pre-registered ASSUMED-parameter sweeps, separate run)",
    "gate_G2_threshold": "fraction >= 0.25 with CI excluding 0",
}
(OUT / "a3_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
