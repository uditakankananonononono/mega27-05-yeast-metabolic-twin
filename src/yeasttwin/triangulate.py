"""SGD inviable triangulation for the over-rescued essential genes.

Cross-checks each over-rescued gene against curated SGD phenotype records
(data/raw/sgd_phenotype_data.tab.gz) for explicit 'inviable' annotations.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw"
RESULTS = ROOT / "results"

COLS = ['feat_name', 'feat_type', 'gene', 'sgdid', 'ref', 'exp_type',
        'mutant_type', 'allele', 'strain', 'observable', 'qualifier',
        'notes', 'reporter', 'condition']


def main():
    ph = pd.read_csv(DATA / "sgd_phenotype_data.tab.gz", sep='\t',
                     header=None, dtype=str)
    ph.columns = COLS[:len(ph.columns)]
    genes = pd.read_csv(RESULTS / "overrescued_genes.csv")['gene']
    sub = ph[ph['feat_name'].isin(genes)]
    inv = sub[sub['observable'].str.contains('inviable', case=False,
                                             na=False)]
    covered = set(inv['feat_name'])
    out = {'n_overrescued': len(genes),
           'n_inviable_in_sgd': len(covered),
           'missing': [g for g in genes if g not in covered]}
    with open(RESULTS / "sgd_inviable_triangulation.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(out)


if __name__ == "__main__":
    main()
