"""RCSB PDB structural coverage of yeast over-rescue pairs."""
import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parents[1] / "results"


def test_structural_coverage_asymmetry():
    d = json.loads((RESULTS / "pdb_pairs.json").read_text())
    assert d["n_essential"] == 16 and d["n_backup"] == 27
    assert d["essential_with_structure"] == 7
    assert d["backup_with_structure"] == 7
    # 74% of backups carry catalytic annotations with no structure
    assert len(d["backups_without_structure"]) == 20
