"""Fourth species: isozyme over-rescue in M. tuberculosis iNJ661.

Essentiality: Sassetti 2003 TraSH calls as tabulated in Griffin 2011
(PLoS Pathog) Table S2 - 613 essential genes.
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
    model = cobra.io.load_json_model(str(DATA / "iNJ661.json"))
    wt = model.slim_optimize()
    ess_df = pd.read_csv(DATA / "mtb_essentiality_griffin2011_tableS2.tsv",
                         sep="\t")
    ess = ess_df[ess_df.sassetti == "essential"].locus.tolist()
    in_model = [g for g in ess if g in model.genes]
    print(f"iNJ661: {len(model.reactions)} rxns, {len(model.genes)} genes, "
          f"WT {wt:.4f}; Sassetti essentials {len(ess)}, in model {len(in_model)}")
    ko = single_gene_deletion(model, gene_list=in_model, processes=1)
    ko = ko.dropna(subset=["growth"])
    ratios = pd.Series({next(iter(r["ids"])): float(r["growth"]) / wt
                        for _, r in ko.iterrows()})
    out = overrescue_audit(ratios, pd.Series(True, index=in_model),
                           model=model)
    out["species"] = "M. tuberculosis (iNJ661) vs Sassetti2003/Griffin2011"
    out["n_essentials_total"] = len(ess)
    out["n_in_model"] = len(in_model)
    with open(RESULTS / "overrescue_audit_mtb.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps({k: v for k, v in out.items()
                      if k != "overrescued_genes"}, indent=2))
    print("over-rescued:", out["overrescued_genes"])


if __name__ == "__main__":
    main()
