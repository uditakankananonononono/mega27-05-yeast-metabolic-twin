"""Per-study abundance replication: strong in E. coli/Mtb, honest weak in yeast."""
import json
from pathlib import Path

import pandas as pd

RESULTS = Path(__file__).resolve().parents[1] / "results"


def test_ecoli_perstudy_replication():
    c = json.loads((RESULTS / "overrescue_abundance_perstudy_combined.json").read_text())
    ec = c["511145"]
    assert ec["n_testable"] >= 10 and ec["n_significant"] >= 8
    assert ec["fisher_combined_p"] < 1e-6


def test_mtb_perstudy_replication():
    c = json.loads((RESULTS / "overrescue_abundance_perstudy_combined.json").read_text())
    mtb = c["83332"]
    assert mtb["n_testable"] >= 5 and mtb["fisher_combined_p"] < 1e-3


def test_yeast_perstudy_honest_weak():
    """The integrated yeast signal (p=0.038) does NOT replicate per-study.
    This honest negative tempers the yeast mechanism claim; bacteria carry it."""
    c = json.loads((RESULTS / "overrescue_abundance_perstudy_combined.json").read_text())
    y = c["4932"]
    assert y["n_testable"] >= 12
    assert y["fisher_combined_p"] > 0.05  # honest: not robust per-study


def test_perstudy_tables():
    for taxid, n_min in [("4932", 15), ("511145", 15), ("83332", 6), ("99287", 3)]:
        df = pd.read_csv(RESULTS / f"overrescue_abundance_perstudy_{taxid}.csv")
        assert len(df) >= n_min
