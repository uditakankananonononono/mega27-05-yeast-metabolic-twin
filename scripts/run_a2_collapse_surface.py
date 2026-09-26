#!/usr/bin/env python3
"""A2: combined-stress collapse surface on the locked 5,000-env LHS sample.

Sliced execution: --start/--end index into the seeded sample (deterministic,
so slices are reproducible and disjoint). WT only. Writes one CSV per slice;
merge after all slices complete.
"""
import argparse
import csv
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from yeasttwin.model import load_model, model_stats  # noqa: E402
from yeasttwin.climate.apply import prepare  # noqa: E402
from yeasttwin.climate.environment import lhs_sample  # noqa: E402
from yeasttwin.climate.objectives import ethanol_yield, is_collapsed  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("--start", type=int, required=True)
p.add_argument("--end", type=int, required=True)
args = p.parse_args()

envs = lhs_sample()[args.start:args.end]
model = load_model()
stats = model_stats(model)
ref = prepare(model)
ref_ethanol = ethanol_yield(model, envs[0].__class__(), ref).ethanol  # reference env

out = ROOT / "results" / "climate" / f"a2_part_{args.start}_{args.end}.csv"
t0 = time.time()
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["idx", "temperature_c", "ethanol_pct", "osmotic_m",
                "nitrogen_frac", "ethanol", "feasible", "relative_yield",
                "collapsed", "sbml_sha256"])
    for i, env in enumerate(envs):
        r = ethanol_yield(model, env, ref)
        rel = r.ethanol / ref_ethanol if r.feasible and ref_ethanol > 0 else 0.0
        w.writerow([args.start + i, env.temperature_c, env.ethanol_pct,
                    env.osmotic_m, env.nitrogen_frac, r.ethanol, r.feasible,
                    rel, is_collapsed(r.ethanol, r.feasible, ref_ethanol),
                    stats["sbml_sha256"]])
print(f"slice {args.start}-{args.end}: {len(envs)} envs in {time.time()-t0:.1f}s -> {out.name}")
