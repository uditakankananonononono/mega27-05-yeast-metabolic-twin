"""Hermetic tests for the climate package (pre-reg locked behavior).

No network. The in-repo yeast-GEM SBML is the only model input. These tests
pin the locked grid, seeds, parameter-function shapes, constraint mappings
and the reference/collapse definitions so no later edit can silently move
the pre-registered contract.
"""
from __future__ import annotations

import pytest

from yeasttwin.model import load_model
from yeasttwin.climate.environment import (GRID_SEED, SPLIT_SEED, Environment,
                                           REFERENCE, design_test_split,
                                           full_grid, grid_axes, lhs_sample)
from yeasttwin.climate.parameters import (AMMONIUM_EX, GLYCEROL_EX, GROWTH_RXN,
                                          LOCKED)
from yeasttwin.climate.apply import applied, prepare
from yeasttwin.climate.objectives import (COLLAPSE_FRACTION, ethanol_yield,
                                          is_collapsed)


# ---- grid / sampling (locked Section 5) ----

def test_grid_axes_exact():
    ax = grid_axes()
    assert len(ax["temperature_c"]) == 15 and ax["temperature_c"][0] == 28.0 and ax["temperature_c"][-1] == 42.0
    assert len(ax["ethanol_pct"]) == 15 and ax["ethanol_pct"][-1] == 14.0
    assert len(ax["osmotic_m"]) == 16 and abs(ax["osmotic_m"][-1] - 1.5) < 1e-9
    assert len(ax["nitrogen_frac"]) == 20 and abs(ax["nitrogen_frac"][0] - 0.05) < 1e-9
    assert len(full_grid()) == 72000


def test_lhs_reproducible_and_in_range():
    a = lhs_sample()
    b = lhs_sample()
    assert a == b and len(a) == 5000
    c = lhs_sample(seed=1)
    assert c != a
    for env in a[:200]:
        assert 28.0 <= env.temperature_c <= 42.0
        assert 0.0 <= env.ethanol_pct <= 14.0
        assert 0.0 <= env.osmotic_m <= 1.5
        assert 0.05 <= env.nitrogen_frac <= 1.0


def test_split_reproducible_disjoint_complete():
    envs = lhs_sample()
    d1, t1 = design_test_split(envs)
    d2, t2 = design_test_split(envs)
    assert d1 == d2 and t1 == t2
    assert len(d1) == 3500 and len(t1) == 1500
    assert len(set(d1) & set(t1)) == 0
    assert len(set(d1) | set(t1)) == len(set(envs))


# ---- parameter function shapes (locked Section 4) ----

def test_temperature_factor_shape():
    f = LOCKED.temperature.factor
    assert f(30.0) == pytest.approx(1.0)
    assert f(42.0) == 0.0
    assert f(28.0) < 1.0 and f(28.0) > 0.0
    assert 0.0 < f(37.0) < 1.0


def test_ethanol_factor_monotonic():
    f = LOCKED.ethanol.factor
    vals = [f(p) for p in range(0, 15)]
    assert vals[0] == pytest.approx(1.0)
    assert vals[-1] == 0.0
    assert all(vals[i] >= vals[i + 1] for i in range(len(vals) - 1))


# ---- constraint mappings on the real model (locked Sections 4-6) ----

@pytest.fixture(scope="module")
def prepped():
    model = load_model()
    ref = prepare(model)
    return model, ref


def test_reference_env_changes_no_bounds(prepped):
    model, ref = prepped
    before = {r.id: (r.lower_bound, r.upper_bound) for r in model.reactions}
    with applied(model, REFERENCE, ref=ref):
        pass
    after = {r.id: (r.lower_bound, r.upper_bound) for r in model.reactions}
    assert before == after


def test_osmotic_sets_glycerol_lb(prepped):
    model, ref = prepped
    with applied(model, Environment(osmotic_m=1.0), ref=ref):
        assert model.reactions.get_by_id(GLYCEROL_EX).lower_bound == pytest.approx(
            LOCKED.osmotic.k_gly * 1.0)
    assert model.reactions.get_by_id(GLYCEROL_EX).lower_bound == 0.0


def test_nitrogen_scales_ammonium(prepped):
    model, ref = prepped
    with applied(model, Environment(nitrogen_frac=0.2), ref=ref):
        assert model.reactions.get_by_id(AMMONIUM_EX).lower_bound == pytest.approx(-200.0)
    assert model.reactions.get_by_id(AMMONIUM_EX).lower_bound == pytest.approx(-1000.0)


def test_heat_ethanol_cap_growth(prepped):
    model, ref = prepped
    with applied(model, Environment(temperature_c=37.0, ethanol_pct=7.0),
                 ref=ref):
        ub = model.reactions.get_by_id(GROWTH_RXN).upper_bound
        expected = ref.growth_max * LOCKED.temperature.factor(37.0) * LOCKED.ethanol.factor(7.0)
        assert ub == pytest.approx(expected)
        assert ub < ref.growth_max


def test_reference_yield_positive_and_extreme_collapses(prepped):
    model, ref = prepped
    base = ethanol_yield(model, REFERENCE, ref)
    assert base.feasible and base.ethanol > 0
    extreme = ethanol_yield(model, Environment(temperature_c=42.0, ethanol_pct=14.0), ref)
    assert is_collapsed(extreme.ethanol, extreme.feasible, base.ethanol)
    assert not is_collapsed(base.ethanol, base.feasible, base.ethanol)
    assert is_collapsed(0.19 * base.ethanol, True, base.ethanol)
    assert not is_collapsed(COLLAPSE_FRACTION * base.ethanol, True, base.ethanol)  # strict <


def test_heat_alone_reduces_ethanol_robustly(prepped):
    # 37 C: f_T ~ 0.55 -> glucose uptake roughly halves -> ethanol falls well
    # below reference; margin is wide so solver noise cannot flip it.
    model, ref = prepped
    base = ethanol_yield(model, REFERENCE, ref)
    mid = ethanol_yield(model, Environment(temperature_c=37.0), ref)
    assert mid.feasible and 0 < mid.ethanol < 0.9 * base.ethanol


def test_ethanol_axis_caps_fermentation(prepped):
    # 7% v/v: f_E = (1 - 7/14)^1.5 ~ 0.354 -> ethanol ub follows f_E.
    model, ref = prepped
    base = ethanol_yield(model, REFERENCE, ref)
    r = ethanol_yield(model, Environment(ethanol_pct=7.0), ref)
    assert r.feasible and r.ethanol <= 0.36 * base.ethanol * 1.0001


def test_midstress_ngam_rises(prepped):
    from yeasttwin.climate.parameters import NGAM_RXN, NGAM_VALUE
    model, ref = prepped
    with applied(model, Environment(temperature_c=37.0), ref=ref):
        assert model.reactions.get_by_id(NGAM_RXN).lower_bound > NGAM_VALUE


def test_lethal_environment_unviable(prepped):
    model, ref = prepped
    r = ethanol_yield(model, Environment(temperature_c=42.0), ref)
    assert not r.feasible and r.ethanol == 0.0
