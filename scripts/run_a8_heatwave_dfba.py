#!/usr/bin/env python3
"""A8: minimal dynamic heat-wave arm (pre-reg Amendment 3, item 3).

Static-optimization dynamic FBA on the locked medium. TRAJECTORY GRID IS
LOCKED IN THIS FILE BEFORE ANY RUN (Amendment 3: "trajectory grid locked in
code before the run"). No new parameters beyond the locked stress layers;
the only declared constant is the ethanol unit conversion below.

Locked trajectory (16 h, dt = 0.25 h, 64 steps):
  phase 1 baseline      30 C            t in [ 0, 2)
  phase 2 ramp          30 -> 40 C lin  t in [ 2, 6)
  phase 3 hold          40 C            t in [ 6, 8)
  phase 4 peak ramp     40 -> 42 C lin  t in [ 8, 9)
  phase 5 peak hold     42 C            t in [ 9,11)
  phase 6 recovery      42 -> 30 C lin  t in [11,12)
  phase 7 recovered     30 C            t in [12,16]
Control trajectory: constant 30 C for 16 h (same dt), locked here too.

Locked batch setup: working volume 1 L; X0 = 0.1 gDW/L; initial glucose
111 mmol/L (~20 g/L, standard 2% SC); ethanol and glycerol start at 0;
osmotic 0 M; nitrogen 100% (locked base medium). The heat wave is the only
imposed external stress; ethanol stress accumulates ENDOGENOUSLY from the
culture's own secretion - that accumulation is the cumulative exposure.

Declared constant (unit conversion, not a fitted parameter):
  1 % v/v ethanol = 7.894 g/L (ethanol density 0.7894 g/mL at 20 C).

Per step (static optimization): constraints from the locked stress layers
(apply.applied) at the instantaneous (T, accumulated ethanol %v/v); glucose
uptake additionally capped by remaining glucose so the culture cannot
consume more than is present. Objective: maximize growth (dFBA standard);
ethanol/glycerol fluxes recorded at that optimum. If f_T x f_E <
VIABILITY_FRACTION (Amendment 1) or the LP is infeasible, the culture is
stalled: growth and all fluxes 0 for the step (no death modeled - declared
limitation). Euler integration for biomass and metabolites.

Locked comparison vs the static A2 surface (claim scope per Amendment 3:
"cumulative-exposure threshold behavior vs the static A2 surface"):
  - per step, the static collapse test (Section 6) at the instantaneous
    environment: is_collapsed(static max-ethanol flux at same env);
  - collapse onset: first step where the dynamic culture's growth-optimal
    ethanol flux < 20% of the reference max, vs the first step where the
    static surface at the instantaneous env is collapsed;
  - recovery deficit at t = 16 h: 1 - achieved/static-predicted ethanol flux
    at the final instantaneous env (heat-wave) vs control;
  - cumulative exposure: degree-hours above 37 C up to collapse onset.
Outputs: results/climate/a8_heatwave.csv, a8_control.csv, a8_summary.json.
"""
import csv
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from yeasttwin.model import load_model, model_stats  # noqa: E402
from yeasttwin.climate.apply import prepare, applied  # noqa: E402
from yeasttwin.climate.environment import Environment  # noqa: E402
from yeasttwin.climate.objectives import ethanol_yield, is_collapsed  # noqa: E402
from yeasttwin.climate.parameters import (GLUCOSE_EX, ETHANOL_EX, GLYCEROL_EX,  # noqa: E402
                                          GROWTH_RXN, NGAM_RXN, NGAM_VALUE,
                                          BASE_GLUCOSE_LB, VIABILITY_FRACTION,
                                          LOCKED)

# --- locked trajectory grid (do not edit after lock; see module docstring) ---
DT = 0.25            # h
T_END = 16.0         # h
PHASES = [           # (t_start, t_end, T_start, T_end) linear within phase
    (0.0, 2.0, 30.0, 30.0),
    (2.0, 6.0, 30.0, 40.0),
    (6.0, 8.0, 40.0, 40.0),
    (8.0, 9.0, 40.0, 42.0),
    (9.0, 11.0, 42.0, 42.0),
    (11.0, 12.0, 42.0, 30.0),
    (12.0, 16.0 + 1e-12, 30.0, 30.0),
]
CONTROL_T = 30.0
# --- locked batch setup ---
V_L = 1.0
X0 = 0.1             # gDW/L
GLUC0_MMOL = 111.0   # mmol/L (~20 g/L)
ETOH_G_PER_L_PER_PCT = 7.894  # declared unit conversion
MW_GLC = 180.156
MW_ETOH = 46.068
MW_GLY = 92.094


def temperature_at(t: float, control: bool = False) -> float:
    if control:
        return CONTROL_T
    for (t0, t1, T0, T1) in PHASES:
        if t0 <= t <= t1 and (t < t1 or t1 >= T_END):
            frac = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
            return T0 + frac * (T1 - T0)
    return PHASES[-1][3]


def run_trajectory(model, ref, ref_ethanol, control: bool = False):
    rows = []
    X = X0                      # gDW/L
    gluc = GLUC0_MMOL           # mmol/L
    etoh_mmol = 0.0             # mmol/L
    gly_mmol = 0.0              # mmol/L
    t = 0.0
    while t <= T_END + 1e-9:
        T = temperature_at(t, control)
        etoh_g_l = etoh_mmol * MW_ETOH / 1000.0
        etoh_pct = etoh_g_l / ETOH_G_PER_L_PER_PCT
        env = Environment(temperature_c=T, ethanol_pct=etoh_pct,
                          osmotic_m=0.0, nitrogen_frac=1.0)
        f_t = LOCKED.temperature.factor(T)
        f_e = LOCKED.ethanol.factor(etoh_pct)
        stalled = f_t * f_e < VIABILITY_FRACTION or gluc <= 0.0
        mu = u_glc = v_etoh = v_gly = ngam = 0.0
        status = "stalled" if stalled else "ok"
        if not stalled:
            with applied(model, env, LOCKED, ref):
                with model:
                    # cap glucose uptake by what remains in the broth
                    u_cap = gluc / (X * DT)
                    r = model.reactions.get_by_id(GLUCOSE_EX)
                    r.lower_bound = max(r.lower_bound, -u_cap)
                    model.objective = GROWTH_RXN
                    sol = model.optimize()
                    if sol.status == "optimal" and sol.objective_value > 1e-9:
                        mu = float(sol.objective_value)
                        u_glc = -float(sol.fluxes[GLUCOSE_EX])
                        v_etoh = float(sol.fluxes[ETHANOL_EX])
                        v_gly = float(sol.fluxes[GLYCEROL_EX])
                        ngam = float(sol.fluxes[NGAM_RXN])
                    else:
                        status = "infeasible"
        # static A2-surface value at the same instantaneous env
        s = ethanol_yield(model, env, ref)
        static_etoh = s.ethanol
        static_collapsed = is_collapsed(s.ethanol, s.feasible, ref_ethanol)
        dyn_collapsed = (status != "ok") or v_etoh < 0.20 * ref_ethanol
        rows.append(dict(t=t, temperature_c=T, ethanol_pct=etoh_pct,
                         biomass_gdw_l=X, glucose_mmol_l=gluc,
                         ethanol_mmol_l=etoh_mmol, glycerol_mmol_l=gly_mmol,
                         mu=mu, glucose_uptake=u_glc, ethanol_flux=v_etoh,
                         glycerol_flux=v_gly, ngam=ngam, f_t=f_t, f_e=f_e,
                         status=status, static_ethanol=static_etoh,
                         static_collapsed=static_collapsed,
                         dynamic_collapsed=dyn_collapsed))
        # Euler integration
        X += mu * X * DT
        gluc = max(0.0, gluc - u_glc * X * DT)
        etoh_mmol += v_etoh * X * DT
        gly_mmol += v_gly * X * DT
        t += DT
    return rows


def summarize(name, rows, ref_ethanol):
    onset_dyn = next((r["t"] for r in rows if r["dynamic_collapsed"]), None)
    onset_static = next((r["t"] for r in rows if r["static_collapsed"]), None)
    final = rows[-1]
    deg_hours = sum(r["temperature_c"] - 37.0 for r in rows
                    if r["temperature_c"] > 37.0) * DT
    onset_deg_hours = sum(r["temperature_c"] - 37.0 for r in rows
                          if r["temperature_c"] > 37.0
                          and (onset_dyn is None or r["t"] < onset_dyn)) * DT
    rec_deficit = (1.0 - final["ethanol_flux"] / final["static_ethanol"]
                   if final["static_ethanol"] > 0 else None)
    return dict(trajectory=name,
                collapse_onset_dynamic_h=onset_dyn,
                collapse_onset_static_surface_h=onset_static,
                degree_hours_above_37C_total=deg_hours,
                degree_hours_above_37C_before_dynamic_collapse=onset_deg_hours,
                final_biomass_gdw_l=final["biomass_gdw_l"],
                final_ethanol_mmol_l=final["ethanol_mmol_l"],
                final_ethanol_pct=final["ethanol_pct"],
                final_glucose_mmol_l=final["glucose_mmol_l"],
                final_glycerol_mmol_l=final["glycerol_mmol_l"],
                final_ethanol_flux=final["ethanol_flux"],
                final_static_ethanol=final["static_ethanol"],
                recovery_deficit_final=rec_deficit,
                ref_ethanol_max=ref_ethanol)


def main():
    t0 = time.time()
    model = load_model()
    stats = model_stats(model)
    ref = prepare(model)
    ref_ethanol = ethanol_yield(model, Environment(), ref).ethanol
    heat = run_trajectory(model, ref, ref_ethanol, control=False)
    ctrl = run_trajectory(model, ref, ref_ethanol, control=True)
    out = ROOT / "results" / "climate"
    for name, rows in (("a8_heatwave", heat), ("a8_control", ctrl)):
        with open(out / f"{name}.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    summary = dict(
        analysis="A8 minimal dynamic heat-wave dFBA (Amendment 3)",
        locked_grid=dict(dt_h=DT, t_end_h=T_END, phases=PHASES,
                         control_T=CONTROL_T),
        batch=dict(volume_l=V_L, x0_gdw_l=X0, glucose0_mmol_l=GLUC0_MMOL,
                   osmotic_m=0.0, nitrogen_frac=1.0),
        declared_constants=dict(etoh_g_per_l_per_pct=ETOH_G_PER_L_PER_PCT,
                                mw_glc=MW_GLC, mw_etoh=MW_ETOH, mw_gly=MW_GLY),
        sbml_sha256=stats["sbml_sha256"],
        heatwave=summarize("heatwave", heat, ref_ethanol),
        control=summarize("control", ctrl, ref_ethanol),
        runtime_s=round(time.time() - t0, 1))
    with open(out / "a8_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
