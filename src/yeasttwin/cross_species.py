"""Cross-species generalisation test: isozyme over-rescue in E. coli iML1515.

Applies the identical over-rescue audit (Discovery 1) to the E. coli
genome-scale model iML1515 against the Keio essential-gene candidates
(Baba et al. 2006, Supplementary Table 6). Same definitions, same cutoff:
a gene is over-rescued iff experimentally essential, predicted viable
(KO/WT growth > 1e-6), and carrying an isozyme backup in its GPRs.
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


def load_keio_essentials() -> list[str]:
    df = pd.read_excel(DATA / "keio_essential_baba2006_supTable6.xls",
                       header=None, skiprows=4)
    bnums = df[6].dropna().astype(str).str.strip()
    return sorted(b for b in bnums if b.startswith("b") and len(b) == 5)


def main():
    model = cobra.io.load_json_model(str(DATA / "iML1515.json"))
    wt = model.slim_optimize()
    print(f"iML1515: {len(model.reactions)} rxns, {len(model.genes)} genes, "
          f"WT growth {wt:.4f} (medium: default BiGG)")
    keio = load_keio_essentials()
    in_model = [b for b in keio if b in model.genes]
    print(f"Keio essentials: {len(keio)} candidates, {len(in_model)} in model")

    ko = single_gene_deletion(model, gene_list=in_model, processes=1)
    ko = ko.dropna(subset=["growth"])
    ratios = pd.Series({next(iter(r["ids"])): float(r["growth"]) / wt
                        for _, r in ko.iterrows()})

    labels = pd.Series(True, index=in_model)  # all Keio candidates are essentials
    out = overrescue_audit(ratios, labels, model=model)
    out["species"] = "Escherichia coli K-12 (iML1515) vs Keio/Baba2006"
    out["n_keio_candidates"] = len(keio)
    out["n_in_model"] = len(in_model)
    with open(RESULTS / "overrescue_audit_ecoli.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps({k: v for k, v in out.items()
                      if k != "overrescued_genes"}, indent=2))
    print("over-rescued:", out["overrescued_genes"])


if __name__ == "__main__":
    main()
