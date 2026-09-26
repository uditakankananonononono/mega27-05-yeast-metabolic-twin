"""Aggregate Amendment-2 parameter sweeps into a3_parameter_uncertainty.json.

Reads results/climate/sweep_<set>_summary.json for each locked sweep set and
emits the G2-robustness verdict. No model loading; pure JSON aggregation.
"""
import json, glob, os

SETS = ["nominal", "kgly_lo", "kgly_hi", "n_lo", "n_hi", "tmax_lo", "tmax_hi"]
G2_THRESHOLD = 0.25

rows = {}
for s in SETS:
    path = f"results/climate/sweep_{s}_summary.json"
    with open(path) as f:
        rows[s] = json.load(f)

vals = {s: rows[s]["fraction_negative_deviation_in_collapsed"] for s in SETS}
out = {
    "amendment": 2,
    "gate": "G2",
    "gate_threshold": G2_THRESHOLD,
    "a3_primary_5000env": 0.242,
    "a3_primary_ci95": [0.223, 0.261],
    "sweep_envs_per_set": 1000,
    "sets": rows,
    "fraction_negative_deviation_range": [min(vals.values()), max(vals.values())],
    "sets_meeting_G2": [s for s in SETS if vals[s] >= G2_THRESHOLD],
    "sets_not_meeting_G2": [s for s in SETS if vals[s] < G2_THRESHOLD],
    "verdict": ("NOT ROBUST: G2 (>=0.25) holds at nominal, kgly_lo/hi and n_lo/hi "
                "but fails under both temperature-margin perturbations "
                "(tmax_lo, tmax_hi) and on the primary 5000-env A3 run (0.242). "
                "The negative-deviation excess in collapsed environments is real "
                "but its magnitude is sensitive to the thermal-death shape "
                "parameter; it is not a parameter-robust law."),
}
with open("results/climate/a3_parameter_uncertainty.json", "w") as f:
    json.dump(out, f, indent=2)
print(json.dumps({k: out[k] for k in ("fraction_negative_deviation_range", "sets_meeting_G2", "sets_not_meeting_G2")}, indent=2))
