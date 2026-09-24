"""Co-expression (GXA) and STRING channel results - incl. honest negatives."""
import json
from pathlib import Path

import pandas as pd

RESULTS = Path(__file__).resolve().parents[1] / "results"


def test_gxa_coexpression_honest_negative():
    """Yeast backup pairs are NOT less co-expressed than expression-matched
    nulls (median percentile 0.55-0.74): the yeast failure is not a
    co-regulation deficit. Recorded as-is."""
    df = pd.read_csv(RESULTS / "overrescue_coexpression_gxa.csv")
    assert len(df) == 2 and (df.n_pairs == 28).all()
    c = json.loads((RESULTS / "overrescue_coexpression_gxa.json").read_text())
    assert c["fisher_combined_p"] > 0.05


def test_string_yeast_annotation_not_coexpression():
    """Backup pairs are database-channel-linked (median 0.9) but only
    weakly coexpression-linked (median 0.071): annotation transfer."""
    s = json.loads((RESULTS / "string_channels.json").read_text())
    y = s["yeast"]
    assert y["wilcoxon_db_minus_coexpr_p"] < 1e-3
    assert y["median_database"] > 0.8
    assert y["median_coexpression"] < 0.15


def test_string_ecoli_honest_no_signal():
    """E. coli pairs show no database-vs-coexpression contrast (p=0.14):
    the annotation-transfer signature is yeast-specific here."""
    s = json.loads((RESULTS / "string_channels.json").read_text())
    assert s["ecoli"]["wilcoxon_db_minus_coexpr_p"] > 0.05
