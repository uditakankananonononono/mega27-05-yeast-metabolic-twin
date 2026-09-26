"""A4 bottleneck analysis (pre-reg Amendment 2, bounded protocol).

Two locked pieces:
1. Binding-constraint identification at every collapsed environment of the
   Amendment-2 subsample (every 5th env of the 5,000 LHS), objective = max
   ethanol (the A3 collapse objective). Candidate constraints (locked):
   glucose cap r_1714 lb, ethanol cap r_1761 ub, growth cap r_2111 ub,
   NGAM r_4046 (fixed when raised), glycerol drain r_1808 lb, ammonium
   bound r_1654 lb. Binding = flux at bound within BIND_RTOL; duals =
   solver reduced costs where provided.
2. FVA subset, locked here in code BEFORE any FVA outcome is seen:
   subsystems {Glycolysis / gluconeogenesis, Pyruvate metabolism,
   Pentose phosphate pathway, Citrate cycle (TCA cycle), Nitrogen
   metabolism} plus explicit osmoregulatory-glycerol and nitrogen-
   assimilation reactions and the ethanol exchange. FVA at
   fraction_of_optimum = 1.0 (variability AT the optimum, matching the
   binding evidence) on 50 seeded collapsed environments (seed 777).
   FVA support for a constraint = its reaction's flux range width at the
   optimum <= FVA_RTOL of max(1, |flux|). Constraints with lb == ub
   (raised NGAM) are vacuous under FVA and marked so - their evidence is
   the dual value, not FVA.
"""
from __future__ import annotations

import cobra

from .parameters import (AMMONIUM_EX, BASE_AMMONIUM_LB, BASE_GLUCOSE_LB,
                         ETHANOL_EX, GLUCOSE_EX, GLYCEROL_EX, GROWTH_RXN,
                         NGAM_RXN)

BIND_RTOL = 1e-4      # flux within this relative tolerance of a bound = binding
FVA_RTOL = 1e-4       # FVA width <= this fraction of max(1, |flux|) = forced
FVA_SEED = 777
FVA_N_ENVS = 50

FVA_SUBSYSTEMS = (
    "Glycolysis / gluconeogenesis",
    "Pyruvate metabolism",
    "Pentose phosphate pathway",
    "Citrate cycle (TCA cycle)",
    "Nitrogen metabolism",
)
FVA_EXTRA_IDS = (
    # osmoregulatory glycerol module (Hohmann 2002 pathway)
    "r_0164", "r_0487", "r_0488", "r_0489", "r_0490", "r_0491", "r_0492",
    "r_1171", "r_1172", "r_1808", "r_1809",
    # nitrogen assimilation + exchange
    "r_0470", "r_0471", "r_0472", "r_0476", "r_1654",
    # ethanol exchange (objective carrier)
    "r_1761",
    # the constraint-carrier reactions themselves (coverage completeness)
    "r_1714", "r_2111", "r_4046",
)

CONSTRAINTS = (
    # name, reaction id, which bound
    ("glucose_cap", GLUCOSE_EX, "lb"),
    ("ethanol_cap", ETHANOL_EX, "ub"),
    ("growth_cap", GROWTH_RXN, "ub"),
    ("ngam", NGAM_RXN, "fixed"),
    ("glycerol_drain", GLYCEROL_EX, "lb"),
    ("ammonium_bound", AMMONIUM_EX, "lb"),
)


def fva_subset(model: cobra.Model) -> list[str]:
    """Resolve the locked FVA reaction list deterministically."""
    ids = {r.id for r in model.reactions
           if (getattr(r, "subsystem", "") or "") in FVA_SUBSYSTEMS}
    ids.update(FVA_EXTRA_IDS)
    return sorted(ids)


def binding_report(model: cobra.Model) -> dict:
    """Optimize (objective already set) and report binding + duals.

    Returns {name: {"flux","bound","binding","reduced_cost","vacuous_fva"}}.
    Caller must be inside applied(model, env, params, ref).
    """
    sol = model.optimize()
    if sol.status != "optimal":
        return {"_status": sol.status}
    rc = sol.reduced_costs
    out = {"_status": "optimal", "_objective": float(sol.objective_value)}
    for name, rid, kind in CONSTRAINTS:
        r = model.reactions.get_by_id(rid)
        v = float(sol.fluxes[rid])
        if kind == "lb":
            bound = r.lower_bound
            binding = abs(v - bound) <= BIND_RTOL * max(1.0, abs(bound))
            vacuous = r.lower_bound == r.upper_bound
        elif kind == "ub":
            bound = r.upper_bound
            binding = abs(v - bound) <= BIND_RTOL * max(1.0, abs(bound))
            vacuous = r.lower_bound == r.upper_bound
        else:  # fixed
            bound = r.lower_bound
            binding = True
            vacuous = True
        out[name] = {
            "flux": v, "bound": float(bound), "binding": bool(binding),
            "reduced_cost": float(rc[rid]) if rc is not None else None,
            "vacuous_fva": bool(vacuous),
        }
    return out
