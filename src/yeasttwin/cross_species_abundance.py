"""Cross-species abundance mechanism test: E. coli backup expression.

Replicates the PaxDb backup-abundance analysis (yeast) on E. coli:
are the iML1515-trusted backups of the 13 over-rescued Keio essentials
also lowly expressed? Uses PaxDb v4.2 (511145, whole organism).
"""
from __future__ import annotations

import json
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from .paralogs import isozyme_partners

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw"
RESULTS = ROOT / "results"


def load_paxdb_ecoli() -> pd.Series:
    df = pd.read_csv(DATA / "paxdb_511145_whole_organism_integrated.txt",
                     sep="\t", comment="#",
                     names=["internal_id", "string_id", "abundance"])
    df["orf"] = df["string_id"].str.replace("511145.", "", regex=False)
    return df.set_index("orf")["abundance"].astype(float)


def main():
    ppm = load_paxdb_ecoli()
    model = cobra.io.load_json_model(str(DATA / "iML1515.json"))
    genes = json.load(open(RESULTS / "overrescue_audit_ecoli.json"))[
        "overrescued_genes"]
    rows = []
    for gid in genes:
        for p in sorted(isozyme_partners(model, gid)):
            rows.append({"essential": gid, "backup": p,
                         "essential_ppm": ppm.get(gid, np.nan),
                         "backup_ppm": ppm.get(p, np.nan)})
    df = pd.DataFrame(rows).dropna()
    df["log2_ratio_ess_over_backup"] = np.log2(df.essential_ppm / df.backup_ppm)
    df.to_csv(RESULTS / "overrescue_abundance_ecoli.csv", index=False)
    lower = int((df.backup_ppm < df.essential_ppm).sum())
    finite = df[np.isfinite(df.log2_ratio_ess_over_backup)]
    stat = wilcoxon(finite.log2_ratio_ess_over_backup)
    out = {
        "n_pairs": int(len(df)),
        "n_backup_lower": lower,
        "frac_backup_lower": lower / len(df),
        "median_log2_ratio": float(finite.log2_ratio_ess_over_backup.median()),
        "n_pairs_finite": int(len(finite)),
        "wilcoxon_p": float(stat.pvalue),
        "genome_median_ppm": float(ppm.median()),
        "backup_median_ppm": float(df.backup_ppm.median()),
        "essential_median_ppm": float(df.essential_ppm.median()),
    }
    with open(RESULTS / "overrescue_abundance_ecoli.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
