#!/usr/bin/env python3
"""A6 nutrient qualitative checks (design locked in a6/RUNLOG.md pre-run).

E-GEOD-1723 direction check: under N-limited env (nitrogen_frac=0.25),
the limiting nutrient uptake (ammonium r_1654) should saturate its bound
with positive shadow price; non-limiting uptakes (phosphate r_2005,
glucose r_1146, sulfate) should not saturate / have ~zero shadow price.
E-MTAB-8245 spare-capacity check: flux + shadow-price pattern under
N-limitation vs ref, narrative comparison only (no fitting).
Declared scope limitation: the env mapping implements only the nitrogen
lever (AMMONIUM_EX); P/S limitation is not parameterized (parameters.py
has the exchange IDs but apply.py maps only nitrogen_frac) - declared.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "results" / "climate" / "a6"

from yeasttwin.model import load_model  # noqa: E402
from yeasttwin.climate.apply import prepare, applied  # noqa: E402
from yeasttwin.climate.environment import Environment  # noqa: E402
from yeasttwin.climate.parameters import (GROWTH_RXN, LOCKED, AMMONIUM_EX,
                                          PHOSPHATE_EX)  # noqa: E402

GLUCOSE_EX = "r_1146"


def find_sulfate_ex(model):
    for r in model.exchanges:
        if "sulfate" in (r.name or "").lower():
            return r.id
    return None


def snapshot(model, ref, env, tag):
    with applied(model, env, LOCKED, ref):
        model.objective = GROWTH_RXN
        sol = model.optimize()
        out = {"env": tag, "status": sol.status,
               "growth": float(sol.objective_value)}
        sul = find_sulfate_ex(model)
        for label, rid in (("ammonium", AMMONIUM_EX),
                           ("phosphate", PHOSPHATE_EX),
                           ("glucose", GLUCOSE_EX),
                           ("sulfate", sul)):
            if rid is None:
                out[label] = None
                continue
            r = model.reactions.get_by_id(rid)
            f = float(sol.fluxes[rid])
            sp = float(sol.shadow_prices.get(rid, 0.0)) \
                if sol.shadow_prices is not None else None
            out[label] = dict(rxn=rid, flux=f, lb=r.lower_bound,
                              ub=r.upper_bound, shadow_price=sp,
                              at_bound=abs(f - r.lower_bound) < 1e-6
                              or abs(f - r.upper_bound) < 1e-6)
    return out


def main():
    model = load_model()
    ref = prepare(model)
    rows = [snapshot(model, ref, Environment(), "ref"),
            snapshot(model, ref, Environment(nitrogen_frac=0.25),
                     "N_lim_0.25")]
    res = dict(analysis="A6 nutrient qualitative checks (locked design)",
               declared_scope_limitation=(
                   "Only the nitrogen lever is implemented in the env "
                   "mapping; P/S limitation not parameterized."),
               rows=rows)
    fp = OUT / "a6_nutrient.json"
    with open(fp, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
