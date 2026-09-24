"""UniProt curated-record provenance of over-rescue pairs."""
import json
from pathlib import Path

import pandas as pd

RESULTS = Path(__file__).resolve().parents[1] / "results"


def test_entries_fetched():
    y = pd.read_csv(RESULTS / "uniprot_entries_yeast.csv")
    e = pd.read_csv(RESULTS / "uniprot_entries_ecoli.csv")
    assert len(y) == 39 and y.accession.notna().sum() >= 35
    assert len(e) == 31 and e.accession.notna().sum() >= 28


def test_family_similarity_provenance():
    """Nearly all pairs are annotated only by family-level similarity -
    curated-record evidence of annotation transfer."""
    d = json.loads((RESULTS / "uniprot_pairs.json").read_text())
    assert d["yeast"]["both_family_similarity_frac"] > 0.85
    assert d["ecoli"]["both_family_similarity_frac"] > 0.95
    assert d["yeast"]["n_pairs"] == 28 and d["ecoli"]["n_pairs"] == 20
