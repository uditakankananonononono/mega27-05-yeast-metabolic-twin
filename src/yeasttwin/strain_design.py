"""Strain design: rank genetic interventions for succinate overproduction.

Succinate is a well-studied yeast metabolic-engineering target with published
experimental knockout data to validate against (e.g. SDH-complex deletion,
glyoxylate-cycle engagement; Otero et al. 2013, Raab et al. 2010 for S.
cerevisiae succinate production from glucose).

Method: for every single-gene knockout, re-optimize the LP with a mandatory
succinate-exchange lower bound sweep (production-envelope style) and record
the maximum feasible succinate flux at >=10% of wild-type growth. Genes whose
deletion raises attainable succinate flux above the wild-type envelope are
ranked candidates. All LPs run on the in-repo SBML; no live calls.
"""
from __future__ import annotations

import pandas as pd

from .media import complete_y7
from .model import load_model

SUCC_EXCHANGE = "r_2056"  # succinate exchange in yeast-GEM
GROWTH_FLOOR = 0.1


def find_succinate_exchange(model) -> str:
    """Locate the succinate exchange reaction by metabolite name, not id."""
    for m in model.metabolites:
        if m.name.lower().startswith("succinate"):
            for r in m.reactions:
                if r.id.startswith("r_") and len(r.metabolites) == 1:
                    return r.id
    raise ValueError("succinate exchange reaction not found")


def wt_production_envelope(model, points: int = 10) -> pd.DataFrame:
    """Wild-type: max succinate flux at fixed fractions of optimal growth."""
    m = complete_y7(model)
    wt = m.slim_optimize()
    ex = m.reactions.get_by_id(find_succinate_exchange(m))
    rows = []
    for frac in [i / (points - 1) for i in range(points)]:
        with m:
            m.reactions.get_by_id("r_2111").lower_bound = frac * wt
            m.objective = ex
            v = m.slim_optimize(error_value=0.0)
            rows.append({"growth_frac": frac, "succinate_flux": max(v, 0.0)})
    return pd.DataFrame(rows)


def knockout_scan(model, genes: list[str] | None = None,
                  growth_floor: float = GROWTH_FLOOR) -> pd.DataFrame:
    """Max succinate flux per single-gene KO at >= growth_floor * WT growth."""
    m = complete_y7(model)
    wt = m.slim_optimize()
    ex_id = find_succinate_exchange(m)
    genes = genes or [g.id for g in m.genes]
    rows = []
    for gid in genes:
        with m:
            try:
                m.genes.get_by_id(gid).knock_out()
            except KeyError:
                continue
            g = m.slim_optimize(error_value=0.0)
            if g < growth_floor * wt:
                rows.append({"gene": gid, "growth": g,
                             "succinate_flux": 0.0, "viable": False})
                continue
            m.reactions.get_by_id("r_2111").lower_bound = growth_floor * wt
            m.objective = m.reactions.get_by_id(ex_id)
            s = m.slim_optimize(error_value=0.0)
            rows.append({"gene": gid, "growth": g,
                         "succinate_flux": max(s, 0.0), "viable": True})
    return pd.DataFrame(rows).sort_values("succinate_flux", ascending=False)


def growth_coupled_scan(model, genes: list[str] | None = None,
                        coupling: float = 0.9,
                        background: tuple[str, ...] = ()) -> pd.DataFrame:
    """Growth-coupled succinate production per single-gene knockout.

    For each knockout (applied on top of *background* deletions), succinate
    flux is maximised subject to growth >= coupling * the mutant's own
    maximal growth. Wild type scores ~0 under this metric; genes whose
    deletion makes succinate a near-obligatory byproduct of growth rank
    high (the OptKnock notion of growth coupling).
    """
    m = complete_y7(model)
    for bg in background:
        m.genes.get_by_id(bg).knock_out()
    ex_id = find_succinate_exchange(m)
    genes = genes or [g.id for g in m.genes if g.id not in background]
    rows = []
    for gid in genes:
        with m:
            try:
                m.genes.get_by_id(gid).knock_out()
            except KeyError:
                continue
            g = m.slim_optimize(error_value=0.0)
            if g < 1e-6:
                rows.append({"gene": gid, "growth": 0.0,
                             "succinate_coupled": 0.0, "viable": False})
                continue
            m.reactions.get_by_id("r_2111").lower_bound = coupling * g
            m.objective = m.reactions.get_by_id(ex_id)
            s = m.slim_optimize(error_value=0.0)
            rows.append({"gene": gid, "growth": g,
                         "succinate_coupled": max(s, 0.0), "viable": True})
    return pd.DataFrame(rows).sort_values("succinate_coupled", ascending=False)


def envelope_shift_scan(model, genes: list[str] | None = None,
                        growth_frac: float = 0.5,
                        background: tuple[str, ...] = ()) -> pd.DataFrame:
    """Envelope-shift scan at a fixed absolute growth floor.

    For each knockout (on top of *background*), maximize succinate flux
    subject to growth >= growth_frac * WILD-TYPE maximal growth (absolute,
    not relative to the mutant). Mutants too sick to reach the floor are
    scored 0. Ranked deltas vs the wild-type value identify true envelope
    shifters, free of the growth-defect confound of relative coupling.
    """
    m = complete_y7(model)
    wt = m.slim_optimize()
    floor = growth_frac * wt
    # WT reference value at the same floor
    with m:
        m.reactions.get_by_id("r_2111").lower_bound = floor
        m.objective = m.reactions.get_by_id(find_succinate_exchange(m))
        wt_val = max(m.slim_optimize(error_value=0.0), 0.0)
    for bg in background:
        m.genes.get_by_id(bg).knock_out()
    ex_id = find_succinate_exchange(m)
    genes = genes or [g.id for g in m.genes if g.id not in background]
    rows = []
    for gid in genes:
        with m:
            try:
                m.genes.get_by_id(gid).knock_out()
            except KeyError:
                continue
            g = m.slim_optimize(error_value=0.0)
            if g < floor:
                rows.append({"gene": gid, "growth": g, "succinate_at_floor": 0.0,
                             "reaches_floor": False})
                continue
            m.reactions.get_by_id("r_2111").lower_bound = floor
            m.objective = m.reactions.get_by_id(ex_id)
            s = m.slim_optimize(error_value=0.0)
            rows.append({"gene": gid, "growth": g,
                         "succinate_at_floor": max(s, 0.0), "reaches_floor": True})
    df = pd.DataFrame(rows)
    df["shift_vs_wt"] = df["succinate_at_floor"] - wt_val
    return df.sort_values("shift_vs_wt", ascending=False), wt_val


def obligatory_production(model, kos: tuple[str, ...] = (),
                          medium: str = "minimal") -> tuple[float, float]:
    """(min succinate at max growth, max growth) for a KO set.

    Min succinate at (near-)maximal growth is the obligatory production
    level: positive iff succinate excretion is stoichiometrically forced
    by growth (growth-coupled design).
    """
    from .media import complete_y7, minimal_glucose
    m = minimal_glucose(model) if medium == "minimal" else complete_y7(model)
    for k in kos:
        m.genes.get_by_id(k).knock_out()
    g = m.slim_optimize(error_value=0.0)
    if g < 1e-6:
        return 0.0, 0.0
    ex = m.reactions.get_by_id(find_succinate_exchange(m))
    ex.lower_bound, ex.upper_bound = 0.0, 1000.0
    m.reactions.get_by_id("r_2111").lower_bound = 0.999 * g
    m.objective = ex
    m.objective.direction = "min"
    return max(m.slim_optimize(error_value=0.0), 0.0), g
