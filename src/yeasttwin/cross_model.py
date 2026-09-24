"""Cross-reconstruction test: isozyme over-rescue in iMM904 (older yeast GEM).

Same species, same Stanford essentiality labels, different reconstruction:
is over-rescue a yeast-GEM 8.x idiosyncrasy or a generic GPR failure?
"""
from __future__ import annotations

import json
from pathlib import Path

import cobra
import pandas as pd
from cobra.flux_analysis import single_gene_deletion

from .labels import load_labels
from .overrescue import overrescue_audit

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw"
RESULTS = ROOT / "results"


def main():
    model = cobra.io.load_json_model(str(DATA / "iMM904.json"))
    wt = model.slim_optimize()
    labels = load_labels()
    in_model = [g for g in labels.index if g in model.genes]
    ess = [g for g in in_model if labels[g]]
    print(f"iMM904: {len(model.reactions)} rxns, {len(model.genes)} genes, "
          f"WT {wt:.4f}; labels {len(in_model)} in-model, {len(ess)} essential")
    ko = single_gene_deletion(model, gene_list=ess, processes=1)
    ko = ko.dropna(subset=["growth"])
    ratios = pd.Series({next(iter(r["ids"])): float(r["growth"]) / wt
                        for _, r in ko.iterrows()})
    out = overrescue_audit(ratios, pd.Series(True, index=ess), model=model)
    out["model"] = "iMM904 (BiGG) vs Stanford deletion labels"
    with open(RESULTS / "overrescue_audit_imm904.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps({k: v for k, v in out.items()
                      if k != "overrescued_genes"}, indent=2))
    print("over-rescued:", out["overrescued_genes"])


if __name__ == "__main__":
    main()
