"""The four locked objectives and the collapse definition (pre-reg Section 6).

Growth is never an objective; viability is a constraint (Amendment 1:
growth/no-growth boundary f_T x f_E >= VIABILITY_FRACTION, plus NGAM
feasibility). Reference environment: 30 C, 0% EtOH, 0 M osmotic, 100% N.
An environment is collapsed for a genotype when max ethanol flux is below
20% of the reference max, or the run is infeasible/nonviable.
"""
from __future__ import annotations

from dataclasses import dataclass

import cobra

from .apply import Reference, applied
from .environment import Environment, grid_axes
from .parameters import (ETHANOL_EX, VIABILITY_FRACTION, LOCKED,
                         StressParameters)

COLLAPSE_FRACTION = 0.20
TOLERANCE_FRACTION = 0.80


@dataclass
class YieldResult:
    ethanol: float
    feasible: bool


def viable(env: Environment, params: StressParameters = LOCKED) -> bool:
    """Amendment-1 growth/no-growth boundary: f_T x f_E >= VIABILITY_FRACTION."""
    f = params.temperature.factor(env.temperature_c) \
        * params.ethanol.factor(env.ethanol_pct)
    return f >= VIABILITY_FRACTION


def ethanol_yield(model: cobra.Model, env: Environment, ref: Reference,
                  params: StressParameters = LOCKED) -> YieldResult:
    """Max ethanol secretion flux under *env* (model flux units).

    Nonviable environments (Amendment 1) return (0.0, False) without a solve.
    """
    if not viable(env, params):
        return YieldResult(0.0, False)
    with applied(model, env, params, ref):
        with model:
            model.objective = ETHANOL_EX
            sol = model.optimize()
            if sol.status != "optimal":
                return YieldResult(0.0, False)
            return YieldResult(float(max(sol.objective_value, 0.0)), True)


def is_collapsed(ethanol: float, feasible: bool, ref_ethanol: float) -> bool:
    return (not feasible) or ethanol < COLLAPSE_FRACTION * ref_ethanol


def tolerance_scan(model: cobra.Model, axis: str, ref: Reference,
                   params: StressParameters = LOCKED,
                   threshold: float = TOLERANCE_FRACTION) -> float:
    """Extreme grid value on *axis* where ethanol yield stays >= threshold*ref.

    Other axes held at reference. Tolerance-objective primitive for
    temperature, ethanol and nitrogen (Section 6, items 2-4).
    """
    ref_ethanol = ethanol_yield(model, Environment(), ref, params).ethanol
    values = grid_axes()[axis]
    best = float(values[0])
    for v in values:
        env = Environment(**{axis: float(v)})
        r = ethanol_yield(model, env, ref, params)
        if r.feasible and r.ethanol >= threshold * ref_ethanol:
            best = float(v)
        else:
            break
    return best
