"""Locked-gate evaluation: repeated stratified 5-fold gene-held-out CV.

Protocol (locked before looking at ML outcomes):
  * Folds are generated ONCE per repeat from a fixed seed and shared by every
    model, so model comparisons are paired by construction.
  * All preprocessing (feature standardization) is fit on the training fold
    only; classification thresholds are tuned on training-fold scores only.
  * Baseline: the published FBA knockout rule (predict essential when the
    KO/WT growth ratio is low), threshold re-tuned per fold on train genes.
  * Reported metric: MCC per (repeat, fold); delta-MCC of each ML model vs
    the FBA baseline; bootstrap 95% CI of the mean delta over resamples of
    the paired fold-level differences. Gate passes iff CI lower bound > 0.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from . import features as feat_mod
from .labels import load_labels
from .ml import (GCN, SeqCNN, auc_roc, best_threshold, encode_sequences, mcc,
                 normalized_adjacency, train_cnn, train_gcn)
from .model import load_model

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
DATA_RAW = ROOT / "data" / "raw"

N_SPLITS = 5
N_REPEATS = 3
BASE_SEED = 20260924


# ---------------------------------------------------------------- folds
def make_folds(y: np.ndarray, n_splits: int, n_repeats: int, seed: int) -> np.ndarray:
    """fold_id[i, r] = fold of gene i in repeat r; stratified on y."""
    n = len(y)
    folds = np.empty((n, n_repeats), dtype=int)
    for r in range(n_repeats):
        rng = np.random.default_rng(seed + r)
        fid = np.empty(n, dtype=int)
        for cls in np.unique(y):
            idx = np.flatnonzero(y == cls)
            rng.shuffle(idx)
            fid[idx] = np.arange(len(idx)) % n_splits
        folds[:, r] = fid
    return folds


# ---------------------------------------------------------------- data
def load_feature_table() -> pd.DataFrame:
    df = pd.read_csv(RESULTS / "gene_features_partial.csv", index_col=0)
    df.index.name = "gene"
    return df


def build_graph_cache(model=None) -> Path:
    """Cache gene adjacency as npz keyed by model gene order."""
    out = RESULTS / "gene_graph.npz"
    if out.exists():
        return out
    m = model or load_model()
    adj = feat_mod.gene_adjacency(m)
    genes = [g.id for g in m.genes]
    pos = {g: i for i, g in enumerate(genes)}
    src, dst = [], []
    for g, nbrs in adj.items():
        for h in nbrs:
            src.append(pos[g]); dst.append(pos[h])
    np.savez(out, genes=np.array(genes), src=np.array(src), dst=np.array(dst))
    return out


def load_graph(model=None) -> tuple[list[str], torch.Tensor]:
    npz = np.load(build_graph_cache(model))
    genes = [str(g) for g in npz["genes"]]
    src, dst = npz["src"], npz["dst"]
    adj = [[] for _ in genes]
    for a, b in zip(src, dst):
        adj[int(a)].append(int(b))
    return genes, normalized_adjacency(adj, len(genes))


def load_sequences(genes: list[str]) -> list[str]:
    """Map model gene ids to protein sequences via the swissprot alias field.

    The gene_id column carries space-separated aliases (standard names AND
    systematic ORFs, e.g. "RAM2 YKL019W"), so every alias token maps to the
    sequence; first hit wins.
    """
    sw = pd.read_csv(DATA_RAW / "swissprot.tsv", sep="\t")
    by_alias: dict[str, str] = {}
    for _, row in sw.iterrows():
        aliases = str(row.get("gene_id", "")).split()
        seq = str(row["sequence"])
        for a in aliases:
            by_alias.setdefault(a, seq)
    return [by_alias.get(g, "") for g in genes]


# ---------------------------------------------------------------- models
def score_fba(train_feat: pd.DataFrame, test_feat: pd.DataFrame, **_) -> np.ndarray:
    return (1.0 - test_feat["ko_ratio_complete"]).to_numpy()


def score_logreg(train_feat, test_feat, y_train, seed=0) -> np.ndarray:
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    sc = StandardScaler().fit(train_feat)
    clf = LogisticRegression(max_iter=2000, class_weight="balanced")
    clf.fit(sc.transform(train_feat), y_train)
    return clf.predict_proba(sc.transform(test_feat))[:, 1]


@dataclass
class FoldResult:
    model: str
    repeat: int
    fold: int
    mcc: float
    auc: float


def run_cv(n_repeats: int = N_REPEATS, gcn_epochs: int = 300, cnn_epochs: int = 25,
           verbose: bool = True) -> pd.DataFrame:
    labels = load_labels()
    feat = load_feature_table()
    genes = [g for g in feat.index if g in labels.index]
    y = labels.loc[genes].to_numpy(dtype=bool)
    feat = feat.loc[genes]
    feat_cols = [c for c in feat.columns]
    x_raw = feat[feat_cols].to_numpy(dtype=np.float32)

    ggenes, a_hat = load_graph()
    gpos = {g: i for i, g in enumerate(ggenes)}
    order = np.array([gpos[g] for g in genes])  # gene i -> graph node
    x_graph = np.zeros((len(ggenes), x_raw.shape[1]), dtype=np.float32)
    x_graph[order] = x_raw
    y_graph = np.zeros(len(ggenes), dtype=np.float32)
    y_graph[order] = y.astype(np.float32)

    seqs = load_sequences(ggenes)
    tok = encode_sequences(seqs)

    folds = make_folds(y, N_SPLITS, n_repeats, BASE_SEED)
    rows: list[FoldResult] = []

    for r in range(n_repeats):
        for f in range(N_SPLITS):
            te = np.flatnonzero(folds[:, r] == f)
            tr = np.flatnonzero(folds[:, r] != f)
            sc_mean = x_raw[tr].mean(0); sc_std = x_raw[tr].std(0) + 1e-9
            xtr = (x_raw[tr] - sc_mean) / sc_std
            xte = (x_raw[te] - sc_mean) / sc_std
            xg = torch.tensor((x_graph - sc_mean) / sc_std)
            xg_tr = torch.tensor(order[tr]); xg_te = order[te]

            scores: dict[str, np.ndarray] = {}
            # FBA baseline (uses only the KO-ratio column; threshold on train)
            scores["fba_rule"] = (1.0 - (xte[:, 0] * sc_std[0] + sc_mean[0]))
            fba_train_scores = (1.0 - (xtr[:, 0] * sc_std[0] + sc_mean[0]))
            thr_fba = best_threshold(y[tr], fba_train_scores)
            # Logistic regression baseline
            scores["logreg"] = score_logreg(pd.DataFrame(xtr), pd.DataFrame(xte), y[tr])
            # GCN (transductive-masked: labels of test nodes never read)
            s_gcn = train_gcn(xg, a_hat, y_graph, xg_tr, epochs=gcn_epochs,
                              seed=BASE_SEED + r * 100 + f)
            scores["gcn"] = s_gcn[xg_te]
            gcn_train = s_gcn[xg_tr]
            # Sequence CNN
            s_cnn = train_cnn(tok, y_graph, order[tr], epochs=cnn_epochs,
                              seed=BASE_SEED + r * 100 + f)
            scores["cnn"] = s_cnn[xg_te]
            cnn_train = s_cnn[order[tr]]
            # Ensemble: mean of train-standardized GCN+CNN scores
            def z(s, ref):
                return (s - ref.mean()) / (ref.std() + 1e-9)
            scores["gcn_cnn"] = (z(scores["gcn"], gcn_train) + z(scores["cnn"], cnn_train)) / 2

            train_score_map = {
                "fba_rule": (fba_train_scores, thr_fba),
                "logreg": (None, None),
                "gcn": (gcn_train, None),
                "cnn": (cnn_train, None),
                "gcn_cnn": (None, None),
            }
            # thresholds tuned on train predictions only
            thr = {"fba_rule": thr_fba}
            thr["logreg"] = best_threshold(y[tr], scores["logreg"] if False else
                score_logreg(pd.DataFrame(xtr), pd.DataFrame(xtr), y[tr]))
            thr["gcn"] = best_threshold(y[tr], gcn_train)
            thr["cnn"] = best_threshold(y[tr], cnn_train)
            thr["gcn_cnn"] = best_threshold(y[tr], (z(gcn_train, gcn_train) + z(cnn_train, cnn_train)) / 2)

            for name, s in scores.items():
                rows.append(FoldResult(name, r, f, mcc(y[te], s >= thr[name]),
                                       auc_roc(y[te], s)))
            if verbose:
                print(f"repeat {r} fold {f}: " + ", ".join(
                    f"{k} MCC={rows[-len(scores)+i].mcc:.3f}" for i, k in enumerate(scores)),
                    flush=True)
    return pd.DataFrame([asdict(x) for x in rows])


def bootstrap_delta_ci(df: pd.DataFrame, model: str, baseline: str = "fba_rule",
                       n_boot: int = 2000, seed: int = 7) -> tuple[float, float, float]:
    """Paired bootstrap over fold-level MCC deltas (model - baseline)."""
    piv = df.pivot_table(index=["repeat", "fold"], columns="model", values="mcc")
    d = (piv[model] - piv[baseline]).to_numpy()
    rng = np.random.default_rng(seed)
    boots = rng.choice(d, size=(n_boot, len(d)), replace=True).mean(1)
    return float(d.mean()), float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def main():
    RESULTS.mkdir(exist_ok=True)
    df = run_cv()
    df.to_csv(RESULTS / "cv_fold_metrics.csv", index=False)
    summ = df.groupby("model").agg(mcc_mean=("mcc", "mean"), mcc_std=("mcc", "std"),
                                   auc_mean=("auc", "mean"))
    out = {"per_model": summ.to_dict("index"), "delta_vs_fba": {}}
    for m in summ.index:
        if m == "fba_rule":
            continue
        mean, lo, hi = bootstrap_delta_ci(df, m)
        out["delta_vs_fba"][m] = {"mean": mean, "ci95": [lo, hi], "gate_pass": lo > 0}
    (RESULTS / "cv_results.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
