"""OMA ortholog-group recurrence lane: preserves the honest negative.

85 of 107 over-rescued backups map to an OMA group, but no group recurs
in more than 2 species: the over-rescue phenomenon is lineage-specific,
not driven by shared ortholog families across species.
"""
import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results" / "oma_recurrence.json"


def test_oma_mapping_succeeded():
    d = json.loads(RESULTS.read_text())
    assert d["n_genes"] == 107
    assert d["n_mapped"] >= 80  # Cloudflare UA fix: real mapping, not the all-null bug


def test_no_cross_species_ortholog_recurrence():
    d = json.loads(RESULTS.read_text())
    assert d["max_species_per_group"] <= 2  # honest negative: no group spans 3+ species
    assert d["n_groups_recurrent"] <= 3
