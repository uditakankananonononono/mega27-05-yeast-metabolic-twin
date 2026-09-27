#!/usr/bin/env python3
"""A8-EXPLORATORY: VHG-batch variant of the locked A8 heat-wave run.

Per pre-reg Section 12, outcome data generated after the locked A8 run under
a changed setup is EXPLORATORY, never pre-registered. Motivation (declared
after inspecting the locked run): the locked 111 mmol/L glucose batch is
substrate-exhausted at 6.5 h (control) / 8.0 h (heat wave), before the
42 C peak, so the locked run cannot observe peak or recovery behavior.
This variant raises initial glucose to 1,665 mmol/L (~300 g/L), the
very-high-gravity (VHG) bioethanol regime cited in pre-reg Section 5
(Fermentation 7(1):38), keeping substrate available through 16 h.
Trajectory grid, dt, integration, and all stress-layer parameters are
IDENTICAL to the locked run (imported, not re-specified). Osmotic axis
stays 0 for comparability (declared; VHG mash osmolarity is not modeled).
"""
import csv
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import run_a8_heatwave_dfba as a8  # noqa: E402
from yeasttwin.model import load_model, model_stats  # noqa: E402
from yeasttwin.climate.apply import prepare  # noqa: E402
from yeasttwin.climate.environment import Environment  # noqa: E402
from yeasttwin.climate.objectives import ethanol_yield  # noqa: E402

a8.GLUC0_MMOL = 1665.0  # ~300 g/L VHG (300 / 180.156 * 1000)


def main():
    t0 = time.time()
    model = load_model()
    stats = model_stats(model)
    ref = prepare(model)
    ref_ethanol = ethanol_yield(model, Environment(), ref).ethanol
    heat = a8.run_trajectory(model, ref, ref_ethanol, control=False)
    ctrl = a8.run_trajectory(model, ref, ref_ethanol, control=True)
    out = ROOT / "results" / "climate"
    for name, rows in (("a8b_vhg_heatwave", heat), ("a8b_vhg_control", ctrl)):
        with open(out / f"{name}.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    summary = dict(
        analysis="A8 VHG variant - EXPLORATORY (Section 12; locked A8 stands)",
        difference_from_locked="initial glucose 1665 mmol/L (~300 g/L VHG) "
                               "vs locked 111 mmol/L; all else identical",
        glucose0_mmol_l=a8.GLUC0_MMOL,
        sbml_sha256=stats["sbml_sha256"],
        heatwave=a8.summarize("heatwave", heat, ref_ethanol),
        control=a8.summarize("control", ctrl, ref_ethanol),
        runtime_s=round(time.time() - t0, 1))
    with open(out / "a8b_vhg_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
