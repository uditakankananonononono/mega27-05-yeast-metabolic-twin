#!/usr/bin/env python3
"""A6: validation against locked Section-8 datasets (pre-reg Section 8, G1).

ALL operational definitions below are LOCKED in this file and committed
BEFORE any run (Amendment discipline). No parameter is fitted to any
validation outcome.

Locked operational definitions (Section 8 metric is fixed by pre-reg;
these are the pre-committed implementation choices):
- Model sensitivity: for gene g, sens(g) = [gr_KO(stress)/gr_KO(ref)] /
  [gr_WT(stress)/gr_WT(ref)]; SENSITIVE if sens(g) < THETA (0.8).
  gr = growth rate (GROWTH_RXN flux) under applied locked env mapping.
  If gr_KO(ref) <= 1e-6 the gene is unconditionally sick/essential in the
  model -> EXCLUDED from that comparison's universe (recorded). Zero or
  infeasible growth under stress with feasible ref growth -> sens = 0.
- Locked stress exemplars (grid points, declared sub-lethal model-side;
  GSE151784's 24.5% lethal dose exceeds the locked 0-14% grid, declared):
    ethanol : Environment(ethanol_pct=8.0)   (matches 8% screen severity)
    heat    : Environment(temperature_c=39.0)
    sorbitol: Environment(osmotic_m=1.0)
    yp20    : Environment(osmotic_m=1.1)     (20% glucose ~ 1.1 M osmotic)
- Dataset HIT calls (locked thresholds, no peeking at overlap):
    GSE151784: per-ORF mean over 3 reps of log2((post+0.5)/(pre+0.5));
               HIT if mean < -1 (2-fold depletion). Barcode-level rows are
               summed to ORF level (trailing _<digit> barcode suffix
               stripped; -A/-B dubious-ORF suffixes preserved).
    Gibney2013 sd01: HIT if column 'Significant' == 'YES' (authors' call).
    GSE59659: per-ORF per-rep median over that ORF's probes of the series
              matrix log-ratio, then mean over 2 reps; HIT if < -1,
              separately for YP20 and sorbitol. Control probes (no ORF)
              excluded.
- Universe per comparison: model genes intersect dataset ORFs minus
  unconditional-sick exclusions and solver_timeout genes (recorded).
- Metrics (locked Section 8): Fisher exact OR, one-sided 'greater' p;
  BH-FDR 5% across the 4 quantitative comparisons; MCC with 95% CI from
  10,000 bootstrap resamples of the universe, seed 20260927.
- G1 second clause: overlap significant (OR > 1, BH-FDR 5%) in >= 2 of 4
  quantitative comparisons. Reported plainly either way (Section 10).
"""
import gzip
import json
import re
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "results" / "climate" / "a6"
OUT.mkdir(parents=True, exist_ok=True)
RAW = ROOT / "data" / "raw"

from yeasttwin.model import load_model  # noqa: E402
from yeasttwin.climate.apply import prepare, applied  # noqa: E402
from yeasttwin.climate.environment import Environment  # noqa: E402
from yeasttwin.climate.parameters import GROWTH_RXN, LOCKED  # noqa: E402

THETA = 0.8
REGIMES = {
    "ethanol": Environment(ethanol_pct=8.0),
    "heat": Environment(temperature_c=39.0),
    "sorbitol": Environment(osmotic_m=1.0),
    "yp20": Environment(osmotic_m=1.1),
}
GENE_TIMEOUT_S = 20
BOOT_N, BOOT_SEED = 10_000, 20260927

_MOD = None
_REF = None


def _worker_init():
    global _MOD, _REF
    _MOD = load_model()
    _REF = prepare(_MOD)


def _growth(model, ref, env):
    with applied(model, env, LOCKED, ref):
        model.objective = GROWTH_RXN
        sol = model.optimize()
        if sol.status != "optimal" or sol.objective_value is None:
            return 0.0
        return float(sol.objective_value)


def _ko_growth_one(args):
    gene, env = args
    with applied(_MOD, env, LOCKED, _REF):
        with _MOD:
            _MOD.genes.get_by_id(gene).knock_out()
            _MOD.objective = GROWTH_RXN
            sol = _MOD.optimize()
            if sol.status != "optimal" or sol.objective_value is None:
                return 0.0
            return float(sol.objective_value)


def model_sensitivity(genes):
    """Resume-safe: results/climate/a6/growth_<regime>.jsonl per gene."""
    res = {}
    for tag, env in [("ref", Environment())] + list(REGIMES.items()):
        fp = OUT / f"growth_{tag}.jsonl"
        done = {}
        if fp.exists():
            for line in open(fp):
                r = json.loads(line)
                done[r["gene"]] = r
        res[tag] = done
        todo = [g for g in genes if g not in done]
        if not todo:
            continue
        pool = Pool(2, initializer=_worker_init)
        main_mod = load_model()
        main_ref = prepare(main_mod)
        wt = _growth(main_mod, main_ref, env)  # WT growth under env (once)
        del main_mod, main_ref
        def restart():
            nonlocal pool
            pool.terminate(); pool.join()
            pool = Pool(2, initializer=_worker_init)
        with open(fp, "a") as fh:
            pend = []
            for g in todo:
                pend.append((g, pool.apply_async(_ko_growth_one, [(g, env)])))
                if len(pend) == 2 or g is todo[-1]:
                    for gg, ar in pend:
                        try:
                            v = ar.get(timeout=GENE_TIMEOUT_S)
                            row = dict(gene=gg, growth=v)
                        except Exception:
                            row = dict(gene=gg, growth=None,
                                       solver_timeout=True)
                            restart()
                        fh.write(json.dumps(row) + "\n"); fh.flush()
                        res[tag][gg] = row
                    pend = []
        pool.terminate(); pool.join()
        res[tag]["__wt__"] = dict(gene="__wt__", growth=wt)
    return res


def load_gse151784():
    pre, post = {}, {}
    for i in (1, 2, 3):
        for tag, acc in (("pre", f"GSM45909{19+i}"), ("post", f"GSM45909{22+i}")):
            d = pre if tag == "pre" else post
            with gzip.open(RAW / "gse151784" / f"{acc}_{tag}_rep_{i}_counts.txt.gz",
                           "rt") as fh:
                next(fh)
                for line in fh:
                    gene, cnt = line.rstrip("\n").split("\t")
                    gene = re.sub(r"_\d+$", "", gene)
                    d.setdefault(i, {}).setdefault(gene, 0)
                    d[i][gene] += int(cnt)
    hits, vals = set(), {}
    genes = set(pre[1]) | set(post[1])
    for g in genes:
        rs = []
        for i in (1, 2, 3):
            p = pre.get(i, {}).get(g, 0)
            q = post.get(i, {}).get(g, 0)
            rs.append(np.log2((q + 0.5) / (p + 0.5)))
        m = float(np.mean(rs))
        vals[g] = m
        if m < -1.0:
            hits.add(g)
    return hits, vals


def load_gibney():
    hits, vals = set(), {}
    fp = RAW / "gibney2013" / "1318100110_sd01.txt"
    header = None
    for line in open(fp):
        line = line.rstrip("\r\n")
        if not line:
            continue
        f = line.split("\t")
        if header is None:
            header = f
            continue
        orf, sig, dr = f[0], f[3], f[11]
        try:
            vals[orf] = float(dr)
        except ValueError:
            vals[orf] = None
        if sig.strip().upper() == "YES":
            hits.add(orf)
    return hits, vals


def load_gse59659():
    ann = {}
    with gzip.open(RAW / "gse59659" / "GPL9825_old_annotations.txt.gz", "rt") as fh:
        for line in fh:
            if line.startswith(("^", "#")):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "ID":
                continue
            if len(f) > 4 and f[4].strip():
                ann[f[0]] = f[4].strip()
    with gzip.open(RAW / "gse59659" / "GSE59659_series_matrix.txt.gz", "rt") as fh:
        in_tab = False
        header = None
        rows = {"sorb1": {}, "sorb2": {}, "yp1": {}, "yp2": {}}
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_tab = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_tab:
                continue
            f = line.rstrip("\n").replace('"', "").split("\t")
            if header is None:
                header = f
                continue
            pid = f[0]
            orf = ann.get(pid)
            if not orf:
                continue
            for key, col in (("sorb1", 1), ("sorb2", 2), ("yp1", 3), ("yp2", 4)):
                try:
                    rows[key].setdefault(orf, []).append(float(f[col]))
                except ValueError:
                    pass
    def per_orf(keys):
        out = {}
        for orf in set(rows[keys[0]]) | set(rows[keys[1]]):
            vs = [np.median(rows[k][orf]) for k in keys if orf in rows[k]]
            if vs:
                out[orf] = float(np.mean(vs))
        return out
    sorb = per_orf(("sorb1", "sorb2"))
    yp = per_orf(("yp1", "yp2"))
    return ({g for g, v in sorb.items() if v < -1.0}, sorb,
            {g for g, v in yp.items() if v < -1.0}, yp)


def compare(name, model_sens, universe, hits):
    ms = set(model_sens) & universe
    dh = set(hits) & universe
    a = len(ms & dh)
    b = len(ms - dh)
    c = len(dh - ms)
    d = len(universe) - a - b - c
    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    denom = np.sqrt((a+b)*(a+c)*(b+d)*(c+d))
    mcc = ((a*d - b*c) / denom) if denom > 0 else 0.0
    rng = np.random.default_rng(BOOT_SEED)
    uni = np.array(sorted(universe))
    boots = []
    n_ms, n_dh = len(ms), len(dh)
    for _ in range(BOOT_N):
        sm = set(rng.choice(uni, size=n_ms, replace=False))
        sd = set(rng.choice(uni, size=n_dh, replace=False))
        aa = len(sm & sd); bb = n_ms - aa; cc = n_dh - aa
        dd = len(uni) - aa - bb - cc
        den = np.sqrt((aa+bb)*(aa+cc)*(bb+dd)*(cc+dd))
        boots.append((aa*dd - bb*cc) / den if den > 0 else 0.0)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return dict(comparison=name, universe=len(universe),
                model_sensitive=len(ms), dataset_hits=len(dh),
                overlap=a, fisher_or=float(odds) if np.isfinite(odds) else None,
                fisher_p_one_sided=float(p), mcc=float(mcc),
                mcc_ci95=[float(lo), float(hi)])


def main():
    t0 = time.time()
    model = load_model()
    genes = sorted(g.id for g in model.genes)
    del model
    growth = model_sensitivity(genes)
    wt_ref = growth["ref"]["__wt__"]["growth"]
    sens = {}
    excluded = {}
    for tag in REGIMES:
        wt_s = growth[tag]["__wt__"]["growth"]
        wt_ratio = (wt_s / wt_ref) if wt_ref > 0 else None
        ss, ex = set(), {}
        for g in genes:
            r0, r1 = growth["ref"].get(g), growth[tag].get(g)
            if not r0 or not r1 or r0.get("solver_timeout") or r1.get("solver_timeout"):
                ex[g] = "solver_timeout_or_missing"
                continue
            g0, g1 = r0["growth"], r1["growth"]
            if g0 <= 1e-6:
                ex[g] = "unconditional_sick_or_essential"
                continue
            s = (g1 / g0) / wt_ratio if wt_ratio and wt_ratio > 0 else None
            if s is not None and s < THETA:
                ss.add(g)
        sens[tag] = ss
        excluded[tag] = ex
    hits_et, vals_et = load_gse151784()
    hits_hs, vals_hs = load_gibney()
    hits_sb, vals_sb, hits_yp, vals_yp = load_gse59659()
    genes_set = set(genes)
    specs = [
        ("ethanol_gse151784", "ethanol", hits_et, vals_et),
        ("heat_gibney2013", "heat", hits_hs, vals_hs),
        ("sorbitol_gse59659", "sorbitol", hits_sb, vals_sb),
        ("yp20_gse59659", "yp20", hits_yp, vals_yp),
    ]
    comps = []
    for tag, regime, hits, vals in specs:
        uni = (genes_set & set(vals)) - set(excluded[regime])
        comps.append(compare(tag, sens[regime], uni, hits))
    ps = [c["fisher_p_one_sided"] for c in comps]
    order = np.argsort(ps)
    m = len(ps)
    bh = [False] * m
    thr = 0.0
    for rank, idx in enumerate(order, 1):
        if ps[idx] <= 0.05 * rank / m:
            thr = ps[idx]
    for i in range(m):
        bh[i] = ps[i] <= thr
    for c, sig in zip(comps, bh):
        c["bh_fdr5_significant"] = bool(sig)
    n_sig = sum(1 for c in comps
                if c["bh_fdr5_significant"] and (c["fisher_or"] or 0) > 1)
    verdict = dict(
        analysis="A6 Section-8 deletion-screen overlap (locked metrics)",
        theta=THETA, regimes={k: vars(v) for k, v in REGIMES.items()},
        wt_growth_ref=wt_ref,
        comparisons=comps,
        g1_overlap_clause=dict(n_significant_of_4=n_sig, passes=n_sig >= 2),
        runtime_s=round(time.time() - t0, 1))
    with open(OUT / "a6_overlap.json", "w") as fh:
        json.dump(verdict, fh, indent=1)
    print(json.dumps(comps, indent=1))
    print(f"G1 overlap clause: {n_sig}/4 significant; passes={n_sig >= 2}")


if __name__ == "__main__":
    main()
