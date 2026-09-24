"""Enrichment results: g:Profiler (yeast) + KEGG hypergeometric (bacteria)."""
import json
from pathlib import Path

import pandas as pd

RESULTS = Path(__file__).resolve().parents[1] / "results"


def test_gprofiler_yeast():
    df = pd.read_csv(RESULTS / "gprofiler_overrescued_yeast.csv")
    assert len(df) == 8
    assert df.p.min() < 1e-3
    s = json.loads((RESULTS / "gprofiler_summary.json").read_text())
    # over-rescued are functionally impoverished vs caught essentials
    assert s["yeast_caught_contrast"]["n_significant_terms"] > 100
    assert s["yeast_caught_contrast"]["n_genes"] == 65


def test_kegg_enrich_honest_nulls():
    s = json.loads((RESULTS / "kegg_enrich_summary.json").read_text())
    assert s["saureus"]["n_significant"] == 1
    assert any(t["pathway"] == "sau00010" and t["p_bonf"] < 0.05
               for t in s["saureus"]["top"])
    for tag in ["ecoli", "bsubtilis", "mtb", "styphimurium"]:
        assert s[tag]["n_significant"] == 0  # recorded honest nulls
