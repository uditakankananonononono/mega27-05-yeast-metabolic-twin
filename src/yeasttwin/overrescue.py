"""ISOZYME-OVERRESCUE-AUDITOR.

Audits a genome-scale metabolic model for a named failure mode discovered in
this project: genes predicted viable because their only reactions carry an
alternative isozyme clause in the GPR, where the backup does not exist in
vivo ("isozyme over-rescue"). Works on any SBML model plus an essentiality
gold set, so it generalises beyond yeast-GEM.

Usage:
    python -m yeasttwin.overrescue            # audits the in-repo yeast-GEM
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

from .features import gpr_stats
from .labels import load_labels
from .model import load_model

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"


def overrescue_audit(ko_ratios: pd.Series, labels: pd.Series,
                     model=None, viable_cutoff: float = 1e-6) -> dict:
    """Compare isozyme backing between FBA-caught and FBA-missed essentials.

    ko_ratios: KO/WT growth ratio per gene (predicted viable iff > cutoff).
    labels: boolean observed-essential per gene.
    """
    model = model or load_model()
    gpr = gpr_stats(model)
    genes = [g for g in labels.index if g in ko_ratios.index]
    lab = labels.loc[genes]
    pred_viable = ko_ratios.loc[genes] > viable_cutoff
    backed = gpr["min_isozymes"].reindex(genes).fillna(1) > 1

    ess = lab[lab].index
    caught = ess[~pred_viable[ess]]   # predicted inviable, truly essential
    missed = ess[pred_viable[ess]]    # predicted viable, truly essential
    a = int(backed[missed].sum()); b = int((~backed[missed]).sum())
    c = int(backed[caught].sum()); d = int((~backed[caught]).sum())
    odds, p = fisher_exact([[a, b], [c, d]])
    return {
        "n_essential": int(len(ess)),
        "n_missed": int(len(missed)),
        "n_caught": int(len(caught)),
        "missed_isozyme_backed_frac": a / max(len(missed), 1),
        "caught_isozyme_backed_frac": c / max(len(caught), 1),
        "fisher_odds_ratio": float(odds),
        "fisher_p": float(p),
        "overrescued_genes": sorted(backed[missed][backed[missed]].index),
    }


def main():
    ko = pd.read_json(RESULTS / "fba_single_ko.json")["ratio"]
    ko = pd.Series(ko)
    labels = load_labels()
    out = overrescue_audit(ko, labels)
    (RESULTS / "overrescue_audit.json").write_text(json.dumps(out, indent=2))
    pd.Series(out["overrescued_genes"], name="gene").to_csv(
        RESULTS / "overrescued_genes.csv", index=False)
    print(json.dumps({k: v for k, v in out.items() if k != "overrescued_genes"},
                     indent=2))
    print("over-rescued genes:", ", ".join(out["overrescued_genes"]))


if __name__ == "__main__":
    main()
