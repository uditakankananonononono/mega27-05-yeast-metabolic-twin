"""Hermetic test: expression-aware GPR repair benchmark (honest negative)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_gpr_repair_result():
    out = json.loads((ROOT / "results" / "gpr_repair.json").read_text())
    # baseline reproduces the published confusion matrix exactly
    assert (out["baseline"]["tp"], out["baseline"]["fn"],
            out["baseline"]["fp"], out["baseline"]["tn"]) == (65, 94, 15, 933)
    assert abs(out["baseline"]["accuracy"] - 0.9015) < 1e-3
    # the naive blanket repair did NOT help - keep honest
    assert out["delta_mcc"] < 0
    assert out["clauses_dropped"] > 0
