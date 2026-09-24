"""Regression tests locking the strain-design discovery into CI.

Hermetic: LPs on the in-repo SBML only.
"""
import pytest

from yeasttwin.media import complete_y7, minimal_glucose
from yeasttwin.model import load_model
from yeasttwin.strain_design import obligatory_production

RAAB4 = ("YKL148C", "YLL041C", "YNL037C", "YDL066W")  # SDH1 SDH2 IDH1 IDP1


@pytest.fixture(scope="module")
def minimal():
    return minimal_glucose(load_model())


@pytest.fixture(scope="module")
def complete():
    return complete_y7(load_model())


def test_wt_no_obligatory_succinate(minimal):
    s, g = obligatory_production(minimal, (), medium=None)
    assert s == 0.0 and g > 0.8


def test_raab_quadruple_validated_on_minimal(minimal):
    s, g = obligatory_production(minimal, RAAB4, medium=None)
    assert s > 0.05           # obligatory production
    assert g > 0.9 * 0.8386   # growth cost < 10%


def test_quadruple_not_validated_on_complete(complete):
    """Medium dependence: on complete Y7 the quadruple forces nothing."""
    s, _ = obligatory_production(complete, RAAB4, medium=None)
    assert s == 0.0


def test_sdh2_single_beats_quadruple(minimal):
    s_single, g_single = obligatory_production(minimal, ("YLL041C",), medium=None)
    s_quad, g_quad = obligatory_production(minimal, RAAB4, medium=None)
    assert s_single > s_quad
    assert g_single > g_quad


def test_sdh2_ach1_best_design(minimal):
    s, g = obligatory_production(minimal, ("YLL041C", "YBL015W"), medium=None)
    assert s > 0.085          # 0.0892 measured
    assert g > 0.99 * 0.8386  # >=99% of WT growth
