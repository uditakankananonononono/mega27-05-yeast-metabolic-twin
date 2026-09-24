"""Model loading and provenance for the yeast genome-scale metabolic model.

The canonical input is yeast-GEM (SysBioChalmers consensus model of
S. cerevisiae) distributed as SBML, stored in-repo at data/raw/yeast-GEM.xml
so every test and analysis is hermetic and byte-reproducible.
"""
from __future__ import annotations

import hashlib
from functools import lru_cache
from pathlib import Path

import cobra

DATA_RAW = Path(__file__).resolve().parents[2] / "data" / "raw"
DEFAULT_MODEL_PATH = DATA_RAW / "yeast-GEM.xml"


def model_sha256(path: Path = DEFAULT_MODEL_PATH) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(maxsize=1)
def load_model(path: str = str(DEFAULT_MODEL_PATH)) -> cobra.Model:
    return cobra.io.read_sbml_model(path)


def model_stats(model=None) -> dict:
    m = model or load_model()
    return {
        "reactions": len(m.reactions),
        "metabolites": len(m.metabolites),
        "genes": len(m.genes),
        "sbml_sha256": model_sha256(),
    }
