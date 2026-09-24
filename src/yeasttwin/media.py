"""Growth media definitions for yeast-GEM (exchange-reaction bounds).

Kennedy synthetic complete medium ("complete_Y7") mirrors the preset used by
the upstream yeast-GEM essential-gene benchmark (code/modelTests/essentialGenes.m),
so our baseline numbers are directly comparable to the published metrics.
"""
from __future__ import annotations

import cobra

COMPLETE_Y7_CONSTRAINED = (
    "r_1604", "r_1639", "r_1873", "r_1879", "r_1880", "r_1881", "r_1671",
    "r_1883", "r_1757", "r_1891", "r_1889", "r_1810", "r_1993", "r_1893",
    "r_1897", "r_1947", "r_1899", "r_1900", "r_1902", "r_1967",
    "r_1903", "r_1548", "r_1904", "r_2028", "r_2038", "r_1906", "r_2067",
    "r_1911", "r_1912", "r_1913", "r_2090", "r_1914", "r_2106",
)
GLUCOSE_EX = "r_1714"
COMPLETE_Y7_FREE = (
    "r_1672", "r_1654", "r_1992", "r_2005", "r_2060", "r_1861", "r_1832",
    "r_2100", "r_4593", "r_4595", "r_4596", "r_4597", "r_2049", "r_4594",
    "r_4600", "r_2020",
)
# Minimal mineral medium (glucose, ammonium, phosphate, sulfate, O2, ions)
MINIMAL_FREE = COMPLETE_Y7_FREE


def _set_lb(model: cobra.Model, rxn_id: str, lb: float) -> bool:
    if rxn_id in model.reactions:
        model.reactions.get_by_id(rxn_id).lower_bound = lb
        return True
    return False


def complete_y7(model: cobra.Model, glucose: float = 20.0) -> cobra.Model:
    """Return a copy of *model* constrained to Kennedy synthetic complete medium."""
    m = model.copy()
    for rxn in m.exchanges:
        rxn.lower_bound, rxn.upper_bound = 0.0, 1000.0
    for rid in COMPLETE_Y7_CONSTRAINED:
        _set_lb(m, rid, -0.5)
    _set_lb(m, GLUCOSE_EX, -abs(glucose))
    for rid in COMPLETE_Y7_FREE:
        _set_lb(m, rid, -1000.0)
    return m


def minimal_glucose(model: cobra.Model, glucose: float = 10.0, aerobic: bool = True) -> cobra.Model:
    """Minimal mineral medium with glucose as sole carbon source."""
    m = model.copy()
    for rxn in m.exchanges:
        rxn.lower_bound, rxn.upper_bound = 0.0, 1000.0
    for rid in MINIMAL_FREE:
        _set_lb(m, rid, -1000.0)
    _set_lb(m, GLUCOSE_EX, -abs(glucose))
    if not aerobic:
        _set_lb(m, "r_1992", 0.0)  # oxygen exchange
    return m
