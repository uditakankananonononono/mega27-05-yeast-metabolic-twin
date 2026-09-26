#!/usr/bin/env python3
"""Run A1 (single-stress dose-responses) and write results/climate/."""
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from yeasttwin.climate.dose_response import run_a1  # noqa: E402

OUT = ROOT / "results" / "climate"
OUT.mkdir(parents=True, exist_ok=True)

rows = run_a1()
with open(OUT / "a1_dose_response.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
summary = {
    "analysis": "A1 single-stress dose-response (WT)",
    "prereg_sha256": "b2235b022fc90d8791d524e9133507e688480651bd2c7551015f38662816dfc5",
    "n_rows": len(rows),
    "sbml_sha256": rows[0]["sbml_sha256"],
    "reference_relative_yield": 1.0,
}
(OUT / "a1_summary.json").write_text(json.dumps(summary, indent=2))
print(f"A1 done: {len(rows)} rows -> {OUT/'a1_dose_response.csv'}")
