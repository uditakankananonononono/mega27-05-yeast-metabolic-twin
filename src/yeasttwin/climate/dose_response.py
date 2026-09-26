"""A1: single-stress dose-responses for the four locked axes (pre-reg 7.A1).

Wild type only. Each axis is scanned over its full locked grid with the
other three at reference. Output rows feed the independence null used by
A3, so this artifact must exist before any combined-stress run.
"""
from __future__ import annotations

from ..model import load_model, model_stats
from .apply import prepare
from .environment import Environment, grid_axes
from .objectives import ethanol_yield

AXES = ("temperature_c", "ethanol_pct", "osmotic_m", "nitrogen_frac")


def run_a1() -> list[dict]:
    model = load_model()
    stats = model_stats(model)
    ref = prepare(model)
    rows: list[dict] = []
    ref = None
    for axis in AXES:
        for v in grid_axes()[axis]:
            env = Environment(**{axis: float(v)})
            r = ethanol_yield(model, env, ref)
            row = {
                "axis": axis,
                "value": float(v),
                "ethanol": r.ethanol,
                "feasible": r.feasible,
                "sbml_sha256": stats["sbml_sha256"],
            }
            if env == Environment():
                ref = r.ethanol
            rows.append(row)
    if ref is None or ref <= 0:
        raise RuntimeError("reference environment missing or zero-yield in A1")
    for row in rows:
        row["relative_yield"] = row["ethanol"] / ref if row["feasible"] else 0.0
    return rows
