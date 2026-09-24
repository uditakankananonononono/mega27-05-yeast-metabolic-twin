"""Third-species test: isozyme over-rescue in B. subtilis iYO844.

Same audit against Koo 2017 (Cell Systems) essential genes - 257
protein-coding essentials from two genome-scale deletion libraries.
"""
from __future__ import annotations

import json
from pathlib import Path

import cobra
import pandas as pd
from cobra.flux_analysis import single_gene_deletion

from .overrescue import overrescue_audit

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw"
RESULTS = ROOT / "results"


def main():
    model = cobra.io.load_json_model(str(DATA / "iYO844.json"))
    wt = model.slim_optimize()
    ess = pd.read_csv(DATA / "bsub_essential_koo2017.tsv", sep="\t")
    in_model = [b for b in ess.bsu if b in model.genes]
    print(f"iYO844: {len(model.reactions)} rxns, {len(model.genes)} genes, "
          f"WT {wt:.4f}; Koo2017 essentials {len(ess)}, in model {len(in_model)}")
    ko = single_gene_deletion(model, gene_list=in_model, processes=1)
    ko = ko.dropna(subset=["growth"])
    ratios = pd.Series({next(iter(r["ids"])): float(r["growth"]) / wt
                        for _, r in ko.iterrows()})
    labels = pd.Series(True, index=in_model)
    out = overrescue_audit(ratios, labels, model=model)
    out["species"] = "Bacillus subtilis (iYO844) vs Koo2017 deletion libraries"
    out["n_koo_essentials"] = len(ess)
    out["n_in_model"] = len(in_model)
    with open(RESULTS / "overrescue_audit_bsubtilis.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps({k: v for k, v in out.items()
                      if k != "overrescued_genes"}, indent=2))
    print("over-rescued:", out["overrescued_genes"])


if __name__ == "__main__":
    main()
