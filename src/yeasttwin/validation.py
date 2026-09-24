"""External experimental validations of the virtual cell.

Three independent published datasets, all computed live from the in-repo
model (no stored predictions):
  * Tobias 2013 chemostats: measured glucose/O2/NH3 uptakes -> predicted
    growth vs measured growth.
  * Van Hoek 1998 chemostats: measured dilution rate + glucose/O2 uptakes
    -> predicted CO2 production vs measured.
  * Biolog substrate panel: each substrate as sole carbon source ->
    predicted growth/no-growth vs experimental Biolog calls.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .media import minimal_glucose
from .model import load_model

DATA_RAW = Path(__file__).resolve().parents[2] / "data" / "raw"

EX_MAP = {  # nutrient -> exchange reaction in yeast-GEM
    "glucose": "r_1714", "oxygen": "r_1992", "ammonium": "r_1654",
    "carbon dioxide": "r_1672",
}


def _bound(m, ex_id, lb=None, ub=None):
    if ex_id in m.reactions:
        if lb is not None:
            m.reactions.get_by_id(ex_id).lower_bound = lb
        if ub is not None:
            m.reactions.get_by_id(ex_id).upper_bound = ub


def chemostat_growth_r2() -> dict:
    """Tobias 2013: predict growth from measured uptake rates."""
    df = pd.read_csv(DATA_RAW / "chemostatData_Tobias2013.tsv", sep="\t")
    preds, obs = [], []
    for _, r in df.iterrows():
        m = minimal_glucose(load_model())
        _bound(m, EX_MAP["glucose"], lb=-abs(r["GLCxtI"]))
        _bound(m, EX_MAP["oxygen"], lb=-abs(r["O2xtI"]))
        _bound(m, EX_MAP["ammonium"], lb=-abs(r["NH3xtI"]))
        preds.append(max(m.slim_optimize(error_value=0.0), 0.0))
        obs.append(r["experimental growth"])
    preds, obs = np.array(preds), np.array(obs)
    ss_res = float(((preds - obs) ** 2).sum())
    ss_tot = float(((obs - obs.mean()) ** 2).sum())
    return {"dataset": "Tobias2013", "n": len(obs),
            "r2_growth": 1 - ss_res / ss_tot,
            "rmse": float(np.sqrt(ss_res / len(obs)))}


def chemostat_co2_r2() -> dict:
    """Van Hoek 1998: predicted vs measured growth and CO2/glucose ratio.

    The model's maximal growth under the measured uptakes sits ~10% below
    the measured dilution rates, so forcing growth to the dilution rate is
    infeasible; we compare (a) predicted max growth vs Drate and (b) the
    CO2-per-glucose yield at that optimum vs the measured ratio.
    """
    df = pd.read_csv(DATA_RAW / "chemostatData_VanHoek1998.tsv", sep="\t")
    pg, og, pr, oratio = [], [], [], []
    for _, r in df.iterrows():
        m = minimal_glucose(load_model())
        _bound(m, EX_MAP["glucose"], lb=-abs(r["GlucoseUptake"]))
        _bound(m, EX_MAP["oxygen"], lb=-abs(r["O2uptake"]))
        g = max(m.slim_optimize(error_value=0.0), 0.0)
        pg.append(g); og.append(r["Drate"])
        co2 = m.reactions.get_by_id(EX_MAP["carbon dioxide"]).flux
        pr.append(co2 / r["GlucoseUptake"]); oratio.append(r["CO2production"] / r["GlucoseUptake"])
    def r2(a, b):
        a, b = np.array(a), np.array(b)
        ss_res = float(((a - b) ** 2).sum())
        ss_tot = float(((b - b.mean()) ** 2).sum())
        return 1 - ss_res / ss_tot if ss_tot else float("nan")
    return {"dataset": "VanHoek1998", "n": len(og),
            "r2_growth": r2(pg, og), "r2_co2_per_glucose": r2(pr, oratio)}


def biolog_agreement(recompute: bool = True) -> dict:
    """Biolog: growth/no-growth per substrate, model vs experiment."""
    df = pd.read_csv(DATA_RAW / "Biolog_Substrate.tsv", sep="\t")
    agree = int((df["Growth_Biolog"] == df["Growth_Model"]).sum())
    out = {"dataset": "Biolog", "n": len(df),
           "stored_model_accuracy": agree / len(df)}
    if recompute:
        # recompute model calls: sole-carbon minimal medium per substrate
        from .media import MINIMAL_FREE, _set_lb
        calls = []
        base = load_model()
        m = base.copy()
        for rxn in m.exchanges:
            rxn.lower_bound = 0.0
        for rid in MINIMAL_FREE:
            if rid in m.reactions:
                m.reactions.get_by_id(rid).lower_bound = -1000.0
        m.reactions.get_by_id("r_1714").lower_bound = 0.0  # no glucose
        for _, r in df.iterrows():
            nm = str(r["Name_in_Model"]).lower()
            hit = [x for x in m.exchanges
                   if x.name.lower() in (nm, nm + " exchange")]
            if not hit:  # fall back to metabolite-name match
                hit = [x for x in m.exchanges if any(
                    met.name.lower().split("[")[0].strip() == nm
                    for met in x.metabolites)]
            if not hit:
                calls.append(None)
                continue
            with m:
                hit[0].lower_bound = -10.0
                calls.append("G" if m.slim_optimize(error_value=0.0) > 0.01 else "NG")
        ok = [c == b for c, b in zip(calls, df["Growth_Biolog"]) if c]
        out["recomputed_model_accuracy"] = sum(ok) / max(len(ok), 1)
        out["recomputed_n"] = len(ok)
    return out


def main():
    res = [chemostat_growth_r2(), chemostat_co2_r2(), biolog_agreement()]
    df = pd.DataFrame(res)
    df.to_csv(Path(__file__).resolve().parents[2] / "results" /
              "external_validations.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
