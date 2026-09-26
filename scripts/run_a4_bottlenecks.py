#!/usr/bin/env python3
"""A4 bottleneck analysis, Amendment-2 bounded protocol.

--phase binding: all collapsed envs of the locked 1,000-env subsample
(every 5th of the 5,000 LHS), objective = max ethanol. Re-solves each env
and RE-VERIFIES the sweep_nominal.csv collapse flag from the fresh solve
(second-pass live check). Writes a4_binding.csv + a4_binding_summary.json.

--phase fva: FVA on the locked subset (bottlenecks.fva_subset) for 50
seeded collapsed envs (seed 777), fraction_of_optimum=1.0. Checkpointed:
rows append to a4_fva.csv; --budget seconds bounds each invocation.
Writes a4_fva_summary.json when complete.
"""
import argparse, csv, json, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
from cobra.flux_analysis import flux_variability_analysis  # noqa: E402

from yeasttwin.model import load_model  # noqa: E402
from yeasttwin.climate.apply import prepare  # noqa: E402
from yeasttwin.climate.environment import lhs_sample  # noqa: E402
from yeasttwin.climate.objectives import ethanol_yield, viable  # noqa: E402
from yeasttwin.climate.parameters import ETHANOL_EX, LOCKED  # noqa: E402
from yeasttwin.climate.bottlenecks import (BIND_RTOL, FVA_N_ENVS, FVA_RTOL,  # noqa: E402
                                           FVA_SEED, binding_report,
                                           fva_subset, CONSTRAINTS)

RES = ROOT / "results" / "climate"


def collapsed_subsample_envs():
    """(idx, env) rows flagged collapsed in sweep_nominal.csv (the locked
    Amendment-2 subsample = every 5th env of the 5,000 LHS)."""
    envs = lhs_sample()[::5]  # CSV idx = subsample index, NOT full-LHS index
    rows = []
    with open(RES / "sweep_nominal.csv") as f:
        for r in csv.DictReader(f):
            if r["collapsed"] == "True":
                rows.append((int(r["idx"]), envs[int(r["idx"])]))
    return rows


def phase_binding(model, ref):
    rows = collapsed_subsample_envs()
    recs, mismatches = [], 0
    for idx, env in rows:
        if not viable(env):
            recs.append((idx, "nonviable", "", "", "", ""))
            continue
        from yeasttwin.climate.apply import applied
        with applied(model, env, LOCKED, ref):
            with model:
                model.objective = ETHANOL_EX
                rep = binding_report(model)
        if rep["_status"] != "optimal":
            recs.append((idx, rep["_status"], "", "", "", ""))
            continue
        # second-pass live verification of the collapse flag
        fresh = rep["_objective"] < 0.20 * ref.ethanol_max
        if not fresh:
            mismatches += 1
        bound_set = [n for n, _, _ in CONSTRAINTS
                     if n in rep and rep[n]["binding"] and not rep[n]["vacuous_fva"]]
        rcs = {n: rep[n]["reduced_cost"] for n, _, _ in CONSTRAINTS
               if n in rep and rep[n]["reduced_cost"] is not None}
        recs.append((idx, "optimal", rep["_objective"], ";".join(bound_set),
                     json.dumps(rcs), json.dumps(
                         {n: rep[n]["flux"] for n, _, _ in CONSTRAINTS if n in rep})))
    with open(RES / "a4_binding.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["idx", "status", "objective", "binding_set", "reduced_costs", "fluxes"])
        w.writerows(recs)
    from collections import Counter
    c = Counter(r[3] for r in recs if r[1] == "optimal")
    summary = {
        "n_collapsed_input": len(rows),
        "n_optimal": sum(1 for r in recs if r[1] == "optimal"),
        "collapse_flag_mismatches": mismatches,
        "binding_set_counts": dict(c.most_common()),
        "constraint_counts": dict(Counter(
            n for r in recs if r[1] == "optimal" and r[3] for n in r[3].split(";")
        ).most_common()),
    }
    with open(RES / "a4_binding_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("BINDING", json.dumps(summary)[:800])


def phase_fva(model, ref, budget):
    t0 = time.time()
    rows = collapsed_subsample_envs()
    rng = np.random.default_rng(FVA_SEED)
    pick = rng.permutation(len(rows))[:FVA_N_ENVS]
    chosen = sorted(rows[i] for i in pick)
    done = set()
    fva_path = RES / "a4_fva.csv"
    if fva_path.exists():
        with open(fva_path) as f:
            done = {int(r["idx"]) for r in csv.DictReader(f)}
    subset = fva_subset(model)
    from yeasttwin.climate.apply import applied
    with open(fva_path, "a", newline="") as f:
        w = csv.writer(f)
        if not done:
            w.writerow(["idx", "reaction", "minimum", "maximum", "flux_at_opt"])
        for idx, env in chosen:
            if idx in done:
                continue
            if time.time() - t0 > budget:
                print(f"BUDGET-EXIT done={len(done)}")
                return
            if not viable(env):
                w.writerow([idx, "_nonviable", "", "", ""])
                continue
            with applied(model, env, LOCKED, ref):
                with model:
                    model.objective = ETHANOL_EX
                    sol = model.optimize()
                    if sol.status != "optimal":
                        w.writerow([idx, "_infeasible", "", "", ""])
                        continue
                    fva = flux_variability_analysis(
                        model, reaction_list=subset, fraction_of_optimum=1.0)
                    for rid in subset:
                        w.writerow([idx, rid, float(fva.loc[rid, "minimum"]),
                                    float(fva.loc[rid, "maximum"]),
                                    float(sol.fluxes[rid])])
            f.flush()
            done.add(idx)
    # aggregate
    import pandas as pd
    df = pd.read_csv(fva_path)
    n_done = df[df.reaction.isin(subset)].idx.nunique()
    recs = {}
    for name, rid, _ in CONSTRAINTS:
        sub = df[(df.reaction == rid)]
        forced = 0
        tot = 0
        for _, r in sub.iterrows():
            width = abs(r["maximum"] - r["minimum"])
            if width <= FVA_RTOL * max(1.0, abs(r["flux_at_opt"])):
                forced += 1
            tot += 1
        recs[name] = {"reaction": rid, "n_envs": tot, "n_forced_at_optimum": forced,
                      "fraction_forced": (forced / tot) if tot else None}
    summary = {"n_envs_fva": n_done, "subset_size": len(subset),
               "subset_ids": subset, "constraints": recs}
    with open(RES / "a4_fva_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("FVA", json.dumps({k: v for k, v in summary.items() if k != "subset_ids"})[:800])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=["binding", "fva"], required=True)
    ap.add_argument("--budget", type=float, default=45)
    a = ap.parse_args()
    model = load_model()
    ref = prepare(model)
    if a.phase == "binding":
        phase_binding(model, ref)
    else:
        phase_fva(model, ref, a.budget)


if __name__ == "__main__":
    main()
