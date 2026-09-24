"""Essentiality labels from the Stanford deletion-collection classifications.

data/raw/essentialGenes.tsv records, per model gene, the upstream yeast-GEM
benchmark outcome on Kennedy synthetic complete medium. The file's classes
use a viable-positive convention (verified against SGD inviable ORFs and the
FBA knockout ratios in results/gene_features_partial.csv):

    file class  meaning
    TP (933)    predicted viable, observed viable       -> non-essential
    FN (15)     predicted inviable, observed viable     -> non-essential
    TN (65)     predicted inviable, observed inviable   -> ESSENTIAL
    FP (94)     predicted viable, observed inviable     -> ESSENTIAL

The observed ground truth is therefore essential = {TN, FP}: 159 genes,
every one of which appears in data/raw/inviable_orfs.txt (SGD). Under the
essential-positive convention of the published benchmark this is exactly
65TP / 94FN / 15FP / 933TN: accuracy 0.9015, MCC 0.532.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

DATA_RAW = Path(__file__).resolve().parents[2] / "data" / "raw"
ESSENTIAL_CLASSES = ("TN", "FP")  # file uses viable-positive convention


def load_labels(path: Path = DATA_RAW / "essentialGenes.tsv") -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    lab = df.set_index("gene")["classification"].isin(ESSENTIAL_CLASSES)
    lab.name = "essential"
    assert lab.sum() == 159, "expected 159 observed-essential genes"
    return lab
