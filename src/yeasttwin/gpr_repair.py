"""Expression-aware GPR repair: does deleting low-expressed backup clauses
improve essentiality prediction?

Hypothesis from the over-rescue discovery: model GPRs trust backup
isozymes the cell barely expresses. Repair: in every OR-clause GPR, drop
clauses whose most abundant gene sits below the PaxDb genome median
(a priori threshold, no tuning). Genes missing from PaxDb are treated as
unknown and their clauses KEPT (conservative). Then rerun the Kennedy-Y7
single-deletion benchmark and compare against the published baseline.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from .abundance import load_paxdb
from .features import knockout_ratios
from .labels import load_labels
from .media import complete_y7
from .model import load_model

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"


def repair_gprs(model, ppm: pd.Series, threshold: float) -> int:
    """Drop OR-clauses whose most abundant gene is below threshold.
    Returns number of clauses dropped."""
    dropped = 0
    for rxn in model.reactions:
        rule = rxn.gene_reaction_rule.strip()
        if " or " not in rule:
            continue
        clauses = re.split(r"\s+or\s+", rule)
        keep = []
        for c in clauses:
            genes = [g for g in re.findall(r"[A-Za-z0-9_\-]+", c.replace("and", " "))
                     if g in model.genes]
            abund = [ppm.get(g) for g in genes]
            measured = [a for a in abund if a is not None]
            if not measured or max(measured) >= threshold:
                keep.append(c)
            else:
                dropped += 1
        if keep and len(keep) < len(clauses):
            rxn.gene_reaction_rule = " or ".join(
                f"({c})" if " and " in c else c for c in keep)
    return dropped


def benchmark(ratios: pd.Series, labels: pd.Series, cutoff: float = 1e-6) -> dict:
    genes = [g for g in labels.index if g in ratios.index]
    lab = labels.loc[genes].astype(bool)
    pred_inviable = ratios.loc[genes] <= cutoff
    tp = int((lab & pred_inviable).sum()); fn = int((lab & ~pred_inviable).sum())
    fp = int((~lab & pred_inviable).sum()); tn = int((~lab & ~pred_inviable).sum())
    acc = (tp + tn) / len(genes)
    mcc = matthews_corrcoef(lab, pred_inviable)
    auc = roc_auc_score(lab, 1.0 - ratios.loc[genes])
    return {"tp": tp, "fn": fn, "fp": fp, "tn": tn,
            "accuracy": acc, "mcc": mcc, "auc": auc, "n": len(genes)}


def main():
    ppm = load_paxdb()
    threshold = float(ppm.median())
    labels = load_labels()

    base = complete_y7(load_model())
    r0 = knockout_ratios(base)
    b0 = benchmark(r0, labels)
    print("baseline:", {k: round(v, 4) if isinstance(v, float) else v
                         for k, v in b0.items()})

    fixed = complete_y7(load_model())
    n_dropped = repair_gprs(fixed, ppm, threshold)
    r1 = knockout_ratios(fixed)
    b1 = benchmark(r1, labels)
    print(f"dropped {n_dropped} low-expression clauses "
          f"(threshold {threshold:.2f} ppm = genome median)")
    print("repaired:", {k: round(v, 4) if isinstance(v, float) else v
                          for k, v in b1.items()})

    out = {"threshold_ppm": threshold, "clauses_dropped": n_dropped,
           "baseline": b0, "repaired": b1,
           "delta_mcc": b1["mcc"] - b0["mcc"],
           "delta_auc": b1["auc"] - b0["auc"]}
    with open(RESULTS / "gpr_repair.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps({k: out[k] for k in
                      ["threshold_ppm", "clauses_dropped",
                       "delta_mcc", "delta_auc"]}, indent=2))


if __name__ == "__main__":
    main()
