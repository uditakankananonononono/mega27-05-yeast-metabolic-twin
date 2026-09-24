"""Isozyme co-expression test: are over-rescued genes co-expressed with backups?

The per-study abundance analysis left yeast mechanism-weak: backup proteins
are not robustly less abundant study-by-study. A sharper failure mode is
co-regulation: a backup rescues only if it is expressed in the SAME
conditions where the essential gene operates. We test this in the two
baseline yeast RNA-seq experiments in Expression Atlas (E-MTAB-8621 cell
cycle on ethanol, E-MTAB-8626 cell cycle on glucose - the only two yeast
baseline experiments in the Atlas catalog of 4,562), correlating each
over-rescued essential gene with its model backup(s) across
condition-group means (Spearman). Each pair is scored as a percentile
against an expression-matched null: correlations of the same essential
with 200 random genes whose mean expression matches the backup (controls
the low-expression low-correlation confound). Per experiment, pairs'
percentiles are tested below 0.5 (one-sided Wilcoxon); experiments are
combined with Fisher's method (eq. 11).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, rankdata, wilcoxon

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw" / "gxa_yeast"
RESULTS = ROOT / "results"
RNG = np.random.default_rng(5)

ACCS = ["E-MTAB-8621", "E-MTAB-8626"]


def load_experiment(acc: str):
    """group-mean log2(TPM+1) matrix indexed by systematic ORF."""
    df = pd.read_csv(DATA / f"{acc}-tpms.tsv", sep="\t", index_col=0)
    df = df.drop(columns=[c for c in df.columns if "Name" in c])
    mat = df.applymap(lambda s: np.mean([float(x) for x in str(s).split(",")])
                      if pd.notna(s) else np.nan)
    mat = np.log2(mat.astype(float) + 1.0)
    mat = mat.loc[mat.notna().all(axis=1)]
    return mat


def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y)
    rx = rx - rx.mean(); ry = ry - ry.mean()
    d = np.sqrt((rx**2).sum() * (ry**2).sum())
    return float((rx * ry).sum() / d) if d > 0 else 0.0


def experiment_result(acc: str, pairs: pd.DataFrame) -> dict:
    mat = load_experiment(acc)
    genes = mat.index
    pcts, rs = [], []
    for _, r in pairs.iterrows():
        e, b = r["gene"], r["backup"]
        if e not in genes or b not in genes:
            continue
        true_r = spearman(mat.loc[e], mat.loc[b])
        mean_b = mat.loc[b].mean()
        pool = genes[(mat.mean(axis=1) - mean_b).abs() <= 0.5]
        pool = pool[(pool != e) & (pool != b)]
        null = [spearman(mat.loc[e], mat.loc[g])
                for g in RNG.choice(pool, size=min(200, len(pool)),
                                    replace=False)]
        pcts.append(float((np.array(null) < true_r).mean()))
        rs.append(true_r)
    pcts, rs = np.array(pcts), np.array(rs)
    p = float(wilcoxon(pcts - 0.5, alternative="less").pvalue) \
        if len(pcts) >= 8 else np.nan
    return {"accession": acc, "n_pairs": int(len(pcts)),
            "median_pair_spearman": float(np.median(rs)),
            "median_null_percentile": float(np.median(pcts)),
            "frac_below_median_null": float((pcts < 0.5).mean()),
            "wilcoxon_p": p}


def main():
    par = pd.read_csv(RESULTS / "overrescued_paralog_identity.csv")
    pairs = pd.DataFrame(
        [{"gene": r["gene"], "backup": p}
         for _, r in par.iterrows()
         for p in str(r["partners"]).split(";") if p and p != "nan"])
    rows = [experiment_result(a, pairs) for a in ACCS
            if (DATA / f"{a}-tpms.tsv").exists()]
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "overrescue_coexpression_gxa.csv", index=False)
    ps = df.wilcoxon_p.dropna().values
    x2 = -2 * np.sum(np.log(ps))
    out = {"n_experiments": int(len(df)),
           "n_significant": int((ps < 0.05).sum()),
           "median_pair_spearman_overall": float(df.median_pair_spearman.median()),
           "fisher_X2": float(x2), "fisher_df": int(2 * len(ps)),
           "fisher_combined_p": float(chi2.sf(x2, 2 * len(ps)))}
    with open(RESULTS / "overrescue_coexpression_gxa.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(df.to_string(index=False))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
