"""MEMOTE benchmark artifacts for yeast-GEM and iML1515."""
import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parents[1] / "results"


def test_memote_reports():
    assert (RESULTS / "memote_yeastgem.html").stat().st_size > 1e6
    assert (RESULTS / "memote_iml1515_full.html").stat().st_size > 1e6
    s = json.loads((RESULTS / "memote_scores.json").read_text())
    assert 0.4 < s["yeast-GEM"]["memote_total_score"] < 0.6
    assert 0.4 < s["iML1515"]["memote_total_score"] < 0.6
