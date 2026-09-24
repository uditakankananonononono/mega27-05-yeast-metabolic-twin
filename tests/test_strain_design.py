"""Hermetic strain-design tests on the in-repo SBML (LPs only, no network)."""
import pytest

from yeasttwin.model import load_model
from yeasttwin.strain_design import (find_succinate_exchange,
                                     knockout_scan, wt_production_envelope)


@pytest.fixture(scope="module")
def model():
    return load_model()


def test_succinate_exchange_found(model):
    assert find_succinate_exchange(model) == "r_2056"


def test_wt_envelope_tradeoff(model):
    env = wt_production_envelope(model, points=3)
    # classic growth/production tradeoff: max succinate at zero growth
    assert env["succinate_flux"].iloc[0] > env["succinate_flux"].iloc[-1]
    assert env["succinate_flux"].iloc[-1] == pytest.approx(0.0, abs=1e-6)


def test_knockout_scan_small_subset(model):
    df = knockout_scan(model, genes=["YLR174W", "YKL148C", "YOR136W"])
    assert set(df.columns) >= {"gene", "growth", "succinate_flux", "viable"}
    assert len(df) == 3
    assert (df["succinate_flux"] >= 0).all()
