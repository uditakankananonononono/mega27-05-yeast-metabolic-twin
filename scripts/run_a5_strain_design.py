#!/usr/bin/env python3
"""A5: strain-design search (pre-reg Section 7 + Amendments 3-4).

All design decisions are locked in docs/PREREGISTRATION.md Amendment 4 and
in this file BEFORE any run. Staged execution: --stage NAME [--chunk I].
Every stage writes its own artifact under results/climate/a5/ so progress
survives sandbox suspension. WT values are computed once (stage wt) and
cached; every design yield is relative to WT in the same environment.
"""
import argparse
import itertools
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "results" / "climate" / "a5"
OUT.mkdir(parents=True, exist_ok=True)

from yeasttwin.model import load_model, model_stats  # noqa: E402
from yeasttwin.climate.apply import prepare, applied  # noqa: E402
from yeasttwin.climate.environment import (Environment, lhs_sample,  # noqa: E402
                                           design_test_split, grid_axes)
from yeasttwin.climate.objectives import ethanol_yield  # noqa: E402
from yeasttwin.climate.parameters import (ETHANOL_EX, GROWTH_RXN, LOCKED)  # noqa: E402

# --- locked subsample seeds (Amendment 4 item 1) ---
SCREEN_SEED, DESIGN_SEED, TEST_SEED = 557, 558, 559
N_SCREEN, N_CONFIRM = 40, 150
# --- locked scaling list (Amendment 4 item 2, stage 4) ---
SCALING = {"r_0491": "GPD1/GPD2", "r_0489": "GPP1/GPP2", "r_1172": "FPS1",
           "r_2115": "ADH", "r_0959": "PDC", "r_0173": "ALD6",
           "r_1166": "HXT", "r_1115": "MEP", "r_0195": "TPS-complex"}
FOLDS = [0.25, 0.5, 2.0, 4.0]
EVO_SEEDS = [42, 123, 2026]
DESIGN_TIMEOUT_S = 45
EVO_POP, EVO_GEN = 24, 12
BASELINE_SEEDS = 30
BASELINE_PER_SEED = 5

_MOD = None  # worker model
_REF = None


def _worker_init():
    global _MOD, _REF
    _MOD = load_model()
    _REF = prepare(_MOD)


def env_key(e):
    return (round(e.temperature_c, 6), round(e.ethanol_pct, 6),
            round(e.osmotic_m, 6), round(e.nitrogen_frac, 6))


def collapse_map():
    """WT collapse status per LHS env from the merged A2 surface."""
    import csv
    m = {}
    with open(ROOT / "results" / "climate" / "a2_collapse_surface.csv") as fh:
        for r in csv.DictReader(fh):
            k = (round(float(r["temperature_c"]), 6),
                 round(float(r["ethanol_pct"]), 6),
                 round(float(r["osmotic_m"]), 6),
                 round(float(r["nitrogen_frac"]), 6))
            m[k] = r["collapsed"] == "True"
    return m


def make_subsamples():
    envs = lhs_sample()
    design, test = design_test_split(envs)
    cm = collapse_map()
    def stratified(pool, n, seed):
        rng = np.random.default_rng(seed)
        coll = [e for e in pool if cm.get(env_key(e), False)]
        ok = [e for e in pool if not cm.get(env_key(e), False)]
        half = n // 2
        a = rng.choice(len(coll), size=min(half, len(coll)), replace=False)
        b = rng.choice(len(ok), size=min(n - half, len(ok)), replace=False)
        return [coll[i] for i in a] + [ok[i] for i in b]
    subs = {"screen40": stratified(design, N_SCREEN, SCREEN_SEED),
            "design150": stratified(design, N_CONFIRM, DESIGN_SEED),
            "test150": stratified(test, N_CONFIRM, TEST_SEED)}
    with open(OUT / "subsamples.json", "w") as fh:
        json.dump({k: [env_key(e) for e in v] for k, v in subs.items()}, fh)
    return subs


def load_subsamples():
    with open(OUT / "subsamples.json") as fh:
        raw = json.load(fh)
    return {k: [Environment(*t) for t in v] for k, v in raw.items()}


def apply_mods(model, mods):
    """mods: list of ('ko', gene_id) or ('scale', rxn_id, fold)."""
    for mod in mods:
        if mod[0] == "ko":
            try:
                model.genes.get_by_id(mod[1]).knock_out()
            except KeyError:
                pass
        else:
            r = model.reactions.get_by_id(mod[1])
            r.lower_bound = r.lower_bound * mod[2]
            r.upper_bound = r.upper_bound * mod[2]


def eval_design_on_envs(model, ref, mods, envs, wt_vals):
    """Mean relative ethanol yield + new-collapse count vs cached WT."""
    ratios, new_coll = [], 0
    for e in envs:
        wt = wt_vals[env_key(e)]
        with applied(model, e, LOCKED, ref):
            with model:
                apply_mods(model, mods)
                model.objective = ETHANOL_EX
                v = float(max(model.slim_optimize(error_value=0.0), 0.0))
        if wt > 0:
            ratios.append(v / wt)
            if v < 0.2 * wt:
                new_coll += 1
    return (float(np.mean(ratios)) if ratios else 0.0, new_coll)


def _screen_one(args):
    mods, envs, wt = args
    mean_rel, new_coll = eval_design_on_envs(_MOD, _REF, mods, envs, wt)
    return (mods, mean_rel, new_coll)


def tolerances(model, ref, mods, own_ref=None):
    """O2/O3/O4 per Amendment 4 item 3: full-axis scans, no early exit.
    Threshold: 80% of the design's own unstressed reference yield."""
    if own_ref is None:
        own_ref = design_ref_yield(model, ref, mods)
    out = {}
    ax = grid_axes()
    for axis, mode in (("temperature_c", "max"), ("ethanol_pct", "max"),
                       ("nitrogen_frac", "min")):
        passing = []
        for v in ax[axis]:
            env = Environment(**{axis: float(v)})
            with applied(model, env, LOCKED, ref):
                with model:
                    apply_mods(model, mods)
                    model.objective = ETHANOL_EX
                    sol = model.optimize()
                    y = float(max(sol.objective_value, 0.0)) \
                        if sol.status == "optimal" else 0.0
            if own_ref > 0 and y >= 0.8 * own_ref:
                passing.append(float(v))
        out[axis] = (max(passing) if mode == "max" else min(passing)) \
            if passing else None
    return out


def design_ref_yield(model, ref, mods):
    with applied(model, Environment(), LOCKED, ref):
        with model:
            apply_mods(model, mods)
            model.objective = ETHANOL_EX
            sol = model.optimize()
            return float(max(sol.objective_value, 0.0)) \
                if sol.status == "optimal" else 0.0


def confirm(mods_list, envs, wt_vals, tag):
    """design150/test150 confirmation + tolerances + own-reference yields."""
    _worker_init()
    rows = []
    for mods in mods_list:
        mean_rel, new_coll = eval_design_on_envs(_MOD, _REF, mods, envs,
                                                 wt_vals)
        own_ref = design_ref_yield(_MOD, _REF, mods)
        tols = {}
        for axis, mode in (("temperature_c", "max"), ("ethanol_pct", "max"),
                           ("nitrogen_frac", "min")):
            passing = []
            for v in grid_axes()[axis]:
                env = Environment(**{axis: float(v)})
                with applied(_MOD, env, LOCKED, _REF):
                    with _MOD:
                        apply_mods(_MOD, mods)
                        _MOD.objective = ETHANOL_EX
                        sol = _MOD.optimize()
                        y = float(max(sol.objective_value, 0.0)) \
                            if sol.status == "optimal" else 0.0
                if own_ref > 0 and y >= 0.8 * own_ref:
                    passing.append(float(v))
            tols[axis] = (max(passing) if mode == "max" else min(passing)) \
                if passing else None
        rows.append(dict(mods=mods, mean_rel=mean_rel, new_collapses=new_coll,
                         own_ref_yield=own_ref, tolerances=tols))
    with open(OUT / f"confirm_{tag}.json", "w") as fh:
        json.dump(rows, fh, indent=1)
    return rows


def wt_stage():
    t0 = time.time()
    subs = make_subsamples()
    _worker_init()
    wt = {}
    for name, envs in subs.items():
        for e in envs:
            r = ethanol_yield(_MOD, e, _REF)
            wt[str(env_key(e))] = r.ethanol if r.feasible else 0.0
    wt_ref = ethanol_yield(_MOD, Environment(), _REF).ethanol
    wt_tols = tolerances(_MOD, _REF, [], own_ref=wt_ref)
    with open(OUT / "wt_cache.json", "w") as fh:
        json.dump(dict(wt=wt, wt_ref_yield=wt_ref, wt_tolerances=wt_tols,
                       runtime_s=round(time.time() - t0, 1)), fh)
    print("WT stage done", wt_ref, wt_tols, f"{time.time()-t0:.0f}s")


def load_wt():
    with open(OUT / "wt_cache.json") as fh:
        d = json.load(fh)
    wt = {eval(k): v for k, v in d["wt"].items()}
    return wt, d["wt_ref_yield"], d["wt_tolerances"]


def screen(mods_list, envs, wt_vals, tag):
    """Resume-safe: completed designs are kept in screen_<tag>.json and
    skipped on re-run; the file is rewritten after every batch."""
    t0 = time.time()
    fp = OUT / f"screen_{tag}.json"
    done = {}
    if fp.exists():
        for r in json.load(open(fp)):
            done[json.dumps(r["mods"])] = r
    todo = [m for m in mods_list if json.dumps(m) not in done]
    rows = list(done.values())
    # per-design hard timeout: GLPK can hang in native code where Python
    # signals do not reach, so each design runs on a worker whose result is
    # awaited with a timeout; a hung worker is terminated and the design is
    # RECORDED as solver_timeout (never silently skipped - honest register).
    pool = Pool(2, initializer=_worker_init)

    def restart():
        nonlocal pool
        pool.terminate()
        pool.join()
        pool = Pool(2, initializer=_worker_init)

    pending = []
    for m in todo:
        pending.append((m, pool.apply_async(_screen_one,
                                            [(m, envs, wt_vals)])))
        if len(pending) == 2 or m is todo[-1]:
            for mm, ar in pending:
                try:
                    r = ar.get(timeout=DESIGN_TIMEOUT_S)
                    rows.append(dict(mods=r[0], mean_rel=r[1],
                                     new_collapses=r[2]))
                except Exception:
                    rows.append(dict(mods=mm, mean_rel=None,
                                     new_collapses=None,
                                     solver_timeout=True))
                    restart()
            pending = []
            with open(fp, "w") as fh:
                json.dump(rows, fh)
    pool.terminate()
    pool.join()
    rows.sort(key=lambda r: (r["mean_rel"] is None,
                             -(r["mean_rel"] or 0.0),
                             r["new_collapses"] or 0))
    with open(fp, "w") as fh:
        json.dump(rows, fh)
    out = [(r["mods"], r["mean_rel"], r["new_collapses"]) for r in rows]
    n_to = sum(1 for r in rows if r.get("solver_timeout"))
    print(f"screen {tag}: {len(rows)} designs ({len(todo)} new, "
          f"{n_to} solver-timeout) {time.time()-t0:.0f}s; "
          f"top: {out[0][0]} rel={out[0][1]:.4f}")
    return out


def stage_singles(chunk, nchunks):
    subs = load_subsamples()
    wt, _, _ = load_wt()
    genes = [g.id for g in load_model().genes]
    part = [g for i, g in enumerate(genes) if i % nchunks == chunk]
    screen([[("ko", g)] for g in part], subs["screen40"], wt,
           f"singles_{chunk}")


def stage_evo():
    subs = load_subsamples()
    wt, _, _ = load_wt()
    envs = subs["screen40"]
    # genome: top-30 modifications across confirmed stages
    universe = []
    for tag in ("singles_conf", "doubles_conf", "scaling_conf"):
        p = OUT / f"confirm_{tag}.json"
        if p.exists():
            for r in json.load(open(p)):
                for m in r["mods"]:  # r["mods"] is a list of single mods
                    universe.append(tuple(m))
    universe = list(dict.fromkeys(universe))[:30]
    all_rows = []
    for seed in EVO_SEEDS:
        rng = np.random.default_rng(seed)
        pop = []
        while len(pop) < EVO_POP:
            idx = rng.choice(len(universe), size=3, replace=False)
            pop.append(tuple(universe[int(i)] for i in idx))
        for gen in range(EVO_GEN):
            scored = screen([list(p) for p in pop], envs, wt,
                            f"evo_s{seed}_g{gen}")
            # screen returns (mods, rel, nc) sorted desc by rel
            elite = [tuple(s[0]) for s in scored[:EVO_POP // 3]]
            children = []
            while len(children) < EVO_POP - len(elite):
                a = elite[int(rng.integers(len(elite)))]
                b = elite[int(rng.integers(len(elite)))]
                child = tuple(dict.fromkeys(tuple(a) + tuple(b)))[:3]
                if rng.random() < 0.3:
                    child = tuple(dict.fromkeys(
                        child + (universe[rng.integers(len(universe))],)))[:3]
                children.append(child)
            pop = elite + children
        all_rows.extend(scored[:5])
    uniq = {}
    for mods, rel, nc in all_rows:
        uniq[tuple(mods)] = (rel, nc)
    finalists = sorted(uniq.items(), key=lambda kv: -kv[1][0])[:10]
    rows = confirm([list(m) for m, _ in finalists], load_subsamples()
                   ["design150"], wt, "evo_conf")
    return rows


def stage_baseline():
    subs = load_subsamples()
    wt, _, _ = load_wt()
    genes = [g.id for g in load_model().genes]
    scale_mods = [("scale", rid, f) for rid in SCALING for f in FOLDS]
    designs = []
    for seed in range(BASELINE_SEEDS):
        rng = np.random.default_rng(10_000 + seed)
        for _ in range(BASELINE_PER_SEED):
            picks = rng.choice(len(genes), size=2, replace=False)
            if rng.random() < 0.7:
                designs.append([("ko", genes[int(picks[0])]),
                                ("ko", genes[int(picks[1])])])
            else:
                designs.append([("ko", genes[int(picks[0])]),
                                scale_mods[int(rng.integers(len(scale_mods)))]])
    rows = []
    _worker_init()
    for m in designs:
        rel, nc = eval_design_on_envs(_MOD, _REF, m, subs["design150"], wt)
        rows.append(dict(mods=m, mean_rel=rel, new_collapses=nc))
    with open(OUT / "baseline_random.json", "w") as fh:
        json.dump(rows, fh)
    vals = sorted(r["mean_rel"] for r in rows)
    p95 = vals[int(0.95 * (len(vals) - 1))]
    print(f"baseline n={len(vals)} p95={p95:.4f} median={vals[len(vals)//2]:.4f}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--stage", required=True,
                   choices=["wt", "singles", "singles2", "doubles", "doubles2",
                            "scaling", "evo", "baseline", "test"])
    p.add_argument("--chunk", type=int, default=0)
    p.add_argument("--nchunks", type=int, default=1)
    a = p.parse_args()
    if a.stage == "wt":
        wt_stage()
    elif a.stage == "singles":
        stage_singles(a.chunk, a.nchunks)
    elif a.stage == "singles2":
        subs = load_subsamples()
        wt, _, _ = load_wt()
        rows = []
        for fp in sorted(OUT.glob("screen_singles_*.json")):
            rows.extend(json.load(open(fp)))
        n_genes = len(load_model().genes)
        if len(rows) < n_genes:
            raise SystemExit(f"singles incomplete: {len(rows)}/{n_genes}")
        rows.sort(key=lambda r: (r["mean_rel"] is None, -(r["mean_rel"] or 0.0), r["new_collapses"] or 0))
        top = [r["mods"] for r in rows[:60]]
        confirm(top, subs["design150"], wt, "singles_conf")
    elif a.stage == "doubles":
        subs = load_subsamples()
        wt, _, _ = load_wt()
        conf = json.load(open(OUT / "confirm_singles_conf.json"))
        top30 = [r["mods"][0][1] for r in conf[:30]]
        pairs = [[("ko", g1), ("ko", g2)]
                 for g1, g2 in itertools.combinations(top30, 2)]
        part = [p for i, p in enumerate(pairs) if i % a.nchunks == a.chunk]
        screen(part, subs["screen40"], wt, f"doubles_{a.chunk}")
    elif a.stage == "scaling":
        subs = load_subsamples()
        wt, _, _ = load_wt()
        mods = [[("scale", rid, f)] for rid in SCALING for f in FOLDS]
        screen(mods, subs["screen40"], wt, "scaling")
        rows = json.load(open(OUT / "screen_scaling.json"))
        rows.sort(key=lambda r: (r["mean_rel"] is None, -(r["mean_rel"] or 0.0)))
        confirm([r["mods"] for r in rows[:10]], subs["design150"], wt,
                "scaling_conf")
    elif a.stage == "evo":
        stage_evo()
    elif a.stage == "doubles2":
        subs = load_subsamples()
        wt, _, _ = load_wt()
        rows = []
        for c in range(a.nchunks):
            rows.extend(json.load(open(OUT / f"screen_doubles_{c}.json")))
        rows.sort(key=lambda r: (r["mean_rel"] is None, -(r["mean_rel"] or 0.0), r["new_collapses"] or 0))
        confirm([r["mods"] for r in rows[:10]], subs["design150"], wt,
                "doubles_conf")
    elif a.stage == "baseline":
        stage_baseline()
    elif a.stage == "test":
        stage_test()


def stage_test():
    """G4 verdict stage (Amendment 4 items 3-5). Locked rules:
    - candidates: union of confirmed stage finalists;
    - plausibility (a): a KO design is screened out if KO growth < 1% of WT
      growth in its target-regime exemplar (locked: the collapsed screen40
      env where the design shows its largest relative gain);
    - plausibility (b): Section-8 datasets not yet acquired -> PENDING,
      reported per gene, does not block;
    - evaluation: test150 (touched once) + tolerance scans;
    - G4(a): >= WT on all four objectives and > on at least one;
    - G4(b): O1 > 95th percentile of the random baseline.
    """
    subs = load_subsamples()
    wt, wt_ref, wt_tols = load_wt()
    _worker_init()
    cands = {}
    for tag in ("singles_conf", "doubles_conf", "scaling_conf", "evo_conf"):
        fp = OUT / f"confirm_{tag}.json"
        if fp.exists():
            for r in json.load(open(fp)):
                key = json.dumps(r["mods"])
                if key not in cands or r["mean_rel"] > cands[key]["mean_rel"]:
                    r["stage"] = tag
                    cands[key] = r
    baseline = json.load(open(OUT / "baseline_random.json"))
    b95 = sorted(r["mean_rel"] for r in baseline)[
        int(0.95 * (len(baseline) - 1))]
    out_rows = []
    for r in sorted(cands.values(), key=lambda x: -x["mean_rel"])[:25]:
        mods = r["mods"]
        # target-regime exemplar: collapsed screen40 env of largest rel gain
        screen_envs = subs["screen40"]
        best_env, best_gain = None, -1.0
        for e in screen_envs:
            wtv = wt[env_key(e)]
            if wtv <= 0:
                continue
            mean_rel_single, _ = eval_design_on_envs(_MOD, _REF, mods, [e],
                                                     wt)
            if mean_rel_single > best_gain:
                best_gain, best_env = mean_rel_single, e
        # plausibility (a): essentiality under the exemplar
        essential = None
        if any(m[0] == "ko" for m in mods):
            with applied(_MOD, best_env, LOCKED, _REF):
                with _MOD:
                    _MOD.objective = GROWTH_RXN
                    sol = _MOD.optimize()
                    wt_g = float(sol.objective_value) \
                        if sol.status == "optimal" else 0.0
                with _MOD:
                    apply_mods(_MOD, mods)
                    _MOD.objective = GROWTH_RXN
                    sol = _MOD.optimize()
                    ko_g = float(sol.objective_value) \
                        if sol.status == "optimal" else 0.0
            essential = (wt_g <= 0) or (ko_g < 0.01 * wt_g)
        # test150 evaluation (touched once)
        mean_rel_t, new_coll_t = eval_design_on_envs(
            _MOD, _REF, mods, subs["test150"], wt)
        own_ref = design_ref_yield(_MOD, _REF, mods)
        tols = tolerances(_MOD, _REF, mods, own_ref=own_ref)

        def not_worse(a):
            if tols[a] is None or wt_tols[a] is None:
                return False
            return (tols[a] >= wt_tols[a] if a != "nitrogen_frac"
                    else tols[a] <= wt_tols[a])

        def strictly_better(a):
            return (tols[a] > wt_tols[a] if a != "nitrogen_frac"
                    else tols[a] < wt_tols[a])

        g4a = (mean_rel_t >= 1.0 and all(not_worse(a) for a in tols)
               and (mean_rel_t > 1.0
                    or any(strictly_better(a) for a in tols)))
        g4b = mean_rel_t > b95
        out_rows.append(dict(
            mods=mods, stage=r["stage"], test_mean_rel=mean_rel_t,
            test_new_collapses=new_coll_t, tolerances=tols,
            exemplar_env=env_key(best_env) if best_env else None,
            exemplar_gain=best_gain,
            essential_under_target_regime=essential,
            screened_out=bool(essential),
            filter_b="PENDING (Section-8 datasets not yet acquired)",
            g4a_dominates_wt=g4a, g4b_beats_baseline_p95=g4b,
            baseline_p95=b95))
    verdict = dict(
        analysis="A5 strain design - TEST stage G4 verdict (Amendment 4)",
        wt_tolerances=wt_tols, wt_ref_yield=wt_ref, baseline_p95=b95,
        candidates=out_rows)
    with open(OUT / "test_g4_verdict.json", "w") as fh:
        json.dump(verdict, fh, indent=1)
    n_pass = sum(1 for r in out_rows
                 if r["g4a_dominates_wt"] and r["g4b_beats_baseline_p95"]
                 and not r["screened_out"])
    print(f"G4: {n_pass} of {len(out_rows)} candidates pass all of "
          f"G4(a)+G4(b)+plausibility; baseline p95={b95:.4f}")


if __name__ == "__main__":
    main()
