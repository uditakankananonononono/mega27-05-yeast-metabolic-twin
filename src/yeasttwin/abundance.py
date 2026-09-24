"""Protein-abundance mechanism test for isozyme over-rescue.

Hypothesis: model-trusted backup isozymes cannot rescue in vivo because
they are barely expressed. Uses PaxDb v4.2 whole-organism integrated
abundances (ppm) for S. cerevisiae. Compares each over-rescued essential
gene's abundance against its model backup partner(s).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw"
RESULTS = ROOT / "results"


def load_paxdb() -> pd.Series:
    df = pd.read_csv(DATA / "paxdb_4932_whole_organism_integrated.txt",
                     sep="\t", comment="#",
                     names=["internal_id", "string_id", "abundance"])
    df["orf"] = df["string_id"].str.replace("4932.", "", regex=False)
    return df.set_index("orf")["abundance"].astype(float)


def main():
    ppm = load_paxdb()
    par = pd.read_csv(RESULTS / "overrescued_paralog_identity.csv")
    rows = []
    for _, r in par.iterrows():
        partners = [p for p in str(r["partners"]).split(";") if p and p != "nan"]
        ess_ppm = ppm.get(r["gene"], np.nan)
        for p in partners:
            rows.append({"essential": r["gene"], "backup": p,
                         "essential_ppm": ess_ppm,
                         "backup_ppm": ppm.get(p, np.nan)})
    df = pd.DataFrame(rows).dropna()
    df["log2_ratio_ess_over_backup"] = np.log2(df.essential_ppm / df.backup_ppm)
    df.to_csv(RESULTS / "overrescue_abundance.csv", index=False)

    lower = int((df.backup_ppm < df.essential_ppm).sum())
    stat = wilcoxon(df.log2_ratio_ess_over_backup)
    out = {
        "n_pairs": int(len(df)),
        "n_backup_lower": lower,
        "frac_backup_lower": lower / len(df),
        "median_log2_ratio": float(df.log2_ratio_ess_over_backup.median()),
        "wilcoxon_p": float(stat.pvalue),
        "genome_median_ppm": float(ppm.median()),
        "backup_median_ppm": float(df.backup_ppm.median()),
        "essential_median_ppm": float(df.essential_ppm.median()),
    }
    with open(RESULTS / "overrescue_abundance.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))
    print(df.nsmallest(8, "backup_ppm").to_string(index=False))


if __name__ == "__main__":
    main()
