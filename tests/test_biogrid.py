"""Hermetic test: curated rescue evidence for backup pairs (BioGRID 5.0.261)."""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_backup_rescue_pairs():
    df = pd.read_csv(ROOT / "results" / "biogrid_backup_rescue_pairs.tsv", sep="\t")
    pairs = {frozenset([r.a_sym, r.b_sym]) for r in df.itertuples()}
    assert frozenset(["RER2", "SRT1"]) in pairs
    assert frozenset(["PET9", "AAC3"]) in pairs
    assert len(pairs) >= 4
    assert set(df.system.str.contains("Rescue")) == {True}
