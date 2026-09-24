"""Per-study abundance replication of the over-rescue mechanism.

The integrated PaxDb analysis (abundance.py) showed over-rescued essentials
have higher abundance than the backup isozymes their GPRs trust (yeast
p=0.038, E. coli p=0.007, B. subtilis honest negative p=0.29). Here we
replicate that paired test in EVERY individual proteomics study deposited
in PaxDb v4.2 for each audited organism - 49 per-study datasets, each an
identifier-backed record (taxid-studyID) fetched individually - and combine
per-organism evidence with Fisher's method:

    X^2 = -2 * sum_i ln(p_i),  X^2 ~ chi2(2k)   (Fisher 1932)   [eq. 11]

Organisms: yeast (17 studies), E. coli (18), Mtb (8), B. subtilis (3),
S. Typhimurium (3), H. pylori (3 - audit underpowered, kept for honesty).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
from scipy.stats import chi2, wilcoxon

from .paralogs import isozyme_partners

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw"
PERSTUDY = DATA / "paxdb_perstudy"
RESULTS = ROOT / "results"

# taxid -> (audit json, model json, gff for tag mapping or None)
ORGS = {
    "4932": ("overrescue_audit.json", None, None),          # yeast pairs from CSV
    "511145": ("overrescue_audit_ecoli.json", "iML1515", None),
    "224308": ("overrescue_audit_bsubtilis.json", "iYO844", "gff_bsubtilis_168.gff"),
    "83332": ("overrescue_audit_mtb.json", "iNJ661", "gff_mtb_h37rv.gff"),
    "99287": ("overrescue_audit_styphimurium_deg1011.json", "STM_v1_0",
              "gff_salmonella_lt2.gff"),
}


def gff_tag_map(path: Path):
    """every locus-tag variant -> canonical current locus_tag."""
    canon = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        attrs = dict(re.findall(r"([^=;]+)=([^;]*)", f[8]))
        cur = attrs.get("locus_tag")
        variants = [cur]
        if attrs.get("old_locus_tag"):
            variants += attrs["old_locus_tag"].replace("%2C", ",").split(",")
        for v in variants:
            if v:
                canon[v] = cur
    return canon


def load_study(path: Path, taxid: str, canon) -> pd.Series:
    df = pd.read_csv(path, sep="\t", comment="#", header=None,
                     usecols=[1, 2], names=["string_id", "abundance"])
    df["gene"] = df["string_id"].astype(str).str.replace(f"{taxid}.", "", regex=False)
    if canon:  # translate old STRING locus tags to current ones the model uses
        df["gene"] = df["gene"].map(lambda g: canon.get(g, g))
    return df.set_index("gene")["abundance"].astype(float)


def essential_backup_pairs(taxid: str) -> pd.DataFrame:
    audit, model_id, gff = ORGS[taxid]
    if taxid == "4932":
        par = pd.read_csv(RESULTS / "overrescued_paralog_identity.csv")
        rows = [{"essential": r["gene"], "backup": p}
                for _, r in par.iterrows()
                for p in str(r["partners"]).split(";") if p and p != "nan"]
        return pd.DataFrame(rows)
    model = cobra.io.load_json_model(str(DATA / f"{model_id}.json"))
    genes = json.load(open(RESULTS / audit))["overrescued_genes"]
    rows = [{"essential": g, "backup": b}
            for g in genes for b in sorted(isozyme_partners(model, g))]
    return pd.DataFrame(rows)


def main():
    summary = {}
    for taxid, (audit, model_id, gff) in ORGS.items():
        pairs = essential_backup_pairs(taxid)
        if pairs.empty:
            continue
        canon = gff_tag_map(DATA / gff) if gff else None
        rows = []
        for f in sorted(PERSTUDY.glob(f"{taxid}-*.txt")):
            ppm = load_study(f, taxid, canon)
            d = pairs.copy()
            d["ess"] = d["essential"].map(ppm)
            d["bak"] = d["backup"].map(ppm)
            d = d.dropna()
            d = d[np.isfinite(np.log2(d.ess / d.bak))]
            if len(d) < 6:
                rows.append({"study": f.name, "n_pairs": len(d),
                             "wilcoxon_p": np.nan,
                             "note": "underpowered (<6 pairs)"})
                continue
            stat = wilcoxon(np.log2(d.ess / d.bak))
            rows.append({"study": f.name, "n_pairs": len(d),
                         "median_log2_ratio": float(np.log2(d.ess / d.bak).median()),
                         "frac_backup_lower": float((d.bak < d.ess).mean()),
                         "wilcoxon_p": float(stat.pvalue)})
        df = pd.DataFrame(rows)
        df.to_csv(RESULTS / f"overrescue_abundance_perstudy_{taxid}.csv",
                  index=False)
        ps = df.wilcoxon_p.dropna().values
        x2 = -2 * np.sum(np.log(ps)) if len(ps) else np.nan
        combined_p = float(chi2.sf(x2, 2 * len(ps))) if len(ps) else np.nan
        summary[taxid] = {
            "n_studies": int(len(df)), "n_testable": int(len(ps)),
            "n_significant": int((ps < 0.05).sum()),
            "fisher_X2": float(x2), "fisher_df": int(2 * len(ps)),
            "fisher_combined_p": combined_p,
        }
        print(taxid, json.dumps(summary[taxid]))
    with open(RESULTS / "overrescue_abundance_perstudy_combined.json", "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
