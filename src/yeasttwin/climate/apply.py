"""Apply a combined-stress environment to the model, in place, with restore.

Constraint mappings (pre-registration Section 4 + Amendment 1):
- Heat: glycolytic capacity scales with f_T (temperature-dependent enzyme
  turnover; CTMI form, Rosso 1993; cardinal values anchored to Salvado
  2011/2013) -> glucose uptake lower bound scales by f_T. Growth is
  additionally capped at f_T x f_E x reference growth.
- Ethanol: product inhibition of the fermentation rate (Levenspiel-type,
  anchored to doi:10.1002/bit.260270311) -> ethanol exchange upper bound
  scales by f_E of the reference max.
- Heat x ethanol: maintenance energy NGAM scales inversely with the
  relative growth factor (Amendment 1, ASSUMED form, swept per Section 4).
- Osmotic: obligatory glycerol-production lower bound (osmoregulatory
  carbon drain; Hohmann 2002).
- Nitrogen: ammonium uptake lower bound scales with availability fraction.
- Viability: environments with f_T x f_E < VIABILITY_FRACTION are nonviable
  (Amendment 1 growth/no-growth boundary); declared in objectives, not here.
Setup: prepare(model) sets the locked base medium (complete_y7, Kennedy
synthetic complete matching the upstream essential-gene benchmark;
microaerobic oxygen lb -0.5, documented modeling choice: fully anaerobic
growth is infeasible in this medium, and industrial bioethanol fermentation
is microaerobic) and returns the Reference (unstressed max growth and max
ethanol). applied() mutates only environment-specific bounds and restores
every one on exit, so one working model serves thousands of environments.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass

import cobra

from .. import media
from .environment import Environment
from .parameters import (AMMONIUM_EX, BASE_AMMONIUM_LB, BASE_GLUCOSE_LB,
                         BASE_OXYGEN_LB, ETHANOL_EX, GLUCOSE_EX,
                         GLYCEROL_EX, GROWTH_RXN, NGAM_RXN, NGAM_VALUE,
                         OXYGEN_EX, VIABILITY_FRACTION, LOCKED,
                         StressParameters)


@dataclass(frozen=True)
class Reference:
    """Unstressed optima on the locked base medium (30 C, 0%, 0 M, 100% N)."""
    growth_max: float
    ethanol_max: float


def base_medium(model: cobra.Model) -> None:
    """Constrain *model* in place to the locked fermentation base medium."""
    m = media.complete_y7(model)
    for rxn in m.reactions:
        r = model.reactions.get_by_id(rxn.id)
        r.lower_bound, r.upper_bound = rxn.lower_bound, rxn.upper_bound
    model.reactions.get_by_id(OXYGEN_EX).lower_bound = BASE_OXYGEN_LB


@contextmanager
def applied(model: cobra.Model, env: Environment,
            params: StressParameters = LOCKED,
            ref: Reference | None = None):
    """Yield *model* constrained for *env*; restore all bounds on exit.

    ref: the Reference from prepare(). Pass it for every real run; when
    None (only inside prepare itself) no stress caps are set.
    """
    changed: list[tuple[cobra.Reaction, float, float]] = []

    def set_bounds(rid: str, lb=None, ub=None):
        r = model.reactions.get_by_id(rid)
        changed.append((r, r.lower_bound, r.upper_bound))
        if lb is not None and ub is not None:
            r.bounds = (lb, ub)
        elif lb is not None:
            r.lower_bound = lb
        elif ub is not None:
            r.upper_bound = ub

    try:
        if ref is not None:
            f_t = params.temperature.factor(env.temperature_c)
            f_e = params.ethanol.factor(env.ethanol_pct)
            f = f_t * f_e
            # heat: glycolytic capacity (glucose uptake Vmax proxy)
            if f_t < 1.0:
                set_bounds(GLUCOSE_EX, lb=BASE_GLUCOSE_LB * f_t)
            # ethanol: product inhibition of fermentation rate
            if f_e < 1.0:
                set_bounds(ETHANOL_EX, ub=ref.ethanol_max * f_e)
            # growth cap (heat x ethanol)
            set_bounds(GROWTH_RXN, ub=ref.growth_max * f)
            # Amendment 1: maintenance rises inversely with relative growth;
            # never below the unstressed baseline (clamp at f >= 1)
            ngam = NGAM_VALUE / max(min(f, 1.0), VIABILITY_FRACTION)
            if ngam > NGAM_VALUE:
                set_bounds(NGAM_RXN, lb=ngam, ub=ngam)
        if env.osmotic_m > 0:
            set_bounds(GLYCEROL_EX, lb=params.osmotic.k_gly * env.osmotic_m)
        if env.nitrogen_frac < 1.0:
            set_bounds(AMMONIUM_EX, lb=BASE_AMMONIUM_LB * env.nitrogen_frac)
        yield model
    finally:
        for r, lb, ub in reversed(changed):
            if lb > r.upper_bound:
                r.upper_bound = ub
                r.lower_bound = lb
            else:
                r.lower_bound = lb
                r.upper_bound = ub


def _optimize(model: cobra.Model, rid: str) -> float:
    with model:
        model.objective = rid
        sol = model.optimize()
        if sol.status != "optimal":
            raise RuntimeError(f"reference optimize {rid} infeasible: {sol.status}")
        return float(sol.objective_value)


def prepare(model: cobra.Model) -> Reference:
    """Set the locked base medium in place; return the unstressed Reference."""
    base_medium(model)
    with applied(model, Environment()):
        growth_max = _optimize(model, GROWTH_RXN)
        ethanol_max = _optimize(model, ETHANOL_EX)
        return Reference(growth_max=growth_max, ethanol_max=ethanol_max)
