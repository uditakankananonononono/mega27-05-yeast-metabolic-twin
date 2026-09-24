"""STRING channel decomposition of over-rescue backup pairs.

If backups fail because they are annotation transfers rather than
co-regulated isozymes, over-rescue pairs should score high on STRING's
database/experiment channels (same pathway curation) but LOW on the
co-expression channel (ascore). Tests, per species: (1) paired Wilcoxon
of database-channel minus coexpression-channel score across pairs;
(2) Mann-Whitney of pairs' ascore vs random model-gene pairs' ascore.
Species: yeast (4932) and E. coli (511145). One batched STRING API call
per species (network endpoint, tsv).
"""
from __future__ import annotations

import io
import json
import urllib.parse
import urllib.request
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

from .paralogs import isozyme_partners

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw"
RESULTS = ROOT / "results"
RNG = np.random.default_rng(11)


def string_network(genes: list[str], species: int) -> pd.DataFrame:
    body = urllib.parse.urlencode(
        {"identifiers": "\r".join(genes), "species": species}).encode()
    req = urllib.request.Request(
        "https://string-db.org/api/tsv/network", data=body)
    txt = urllib.request.urlopen(req, timeout=60).read().decode()
    return pd.read_csv(io.StringIO(txt), sep="\t")


def species_result(species: int, pairs: pd.DataFrame,
                   model_genes: list[str]) -> dict:
    genes = sorted(set(pairs.gene) | set(pairs.backup))
    ctrl = RNG.choice([g for g in model_genes if g not in genes],
                      size=min(60, len(model_genes) - len(genes)),
                      replace=False).tolist()
    net = string_network(genes + ctrl, species)
    net["pair"] = net.apply(
        lambda r: tuple(sorted((r.stringId_A.split(".", 1)[1],
                                r.stringId_B.split(".", 1)[1]))), axis=1)
    edge = net.set_index("pair")
    rows, missing = [], 0
    for _, r in pairs.iterrows():
        key = tuple(sorted((r["gene"], r["backup"])))
        if key in edge.index:
            e = edge.loc[[key]].iloc[0]
            rows.append({"species": species, "gene": r["gene"],
                         "backup": r["backup"], "combined": e.score,
                         "coexpression": e.ascore, "experiments": e.escore,
                         "database": e.dscore, "textmining": e.tscore})
        else:
            missing += 1
    df = pd.DataFrame(rows)
    ctrl_pairs = [(tuple(sorted((a, b)))) for i, a in enumerate(ctrl)
                  for b in ctrl[i + 1:]]
    ctrl_a = [edge.loc[[k]].iloc[0].ascore for k in ctrl_pairs
              if k in edge.index]
    out = {"species": species, "n_pairs": int(len(df)),
           "n_missing_edge": missing,
           "median_coexpression": float(df.coexpression.median()),
           "median_database": float(df.database.median()),
           "frac_zero_coexpression": float((df.coexpression == 0).mean()),
           "frac_database_positive": float((df.database > 0).mean()),
           "wilcoxon_db_minus_coexpr_p": float(
               wilcoxon(df.database - df.coexpression,
                        alternative="greater").pvalue),
           "control_n_pairs": len(ctrl_a),
           "control_median_coexpression": float(np.median(ctrl_a)),
           "mannwhitney_pairs_vs_control_p": float(
               mannwhitneyu(df.coexpression, ctrl_a,
                            alternative="less").pvalue)}
    df.to_csv(RESULTS / f"string_channels_{species}.csv", index=False)
    return out


def main():
    par = pd.read_csv(RESULTS / "overrescued_paralog_identity.csv")
    ypairs = pd.DataFrame(
        [{"gene": r["gene"], "backup": p}
         for _, r in par.iterrows()
         for p in str(r["partners"]).split(";") if p and p != "nan"])
    from .model import load_model
    ymodel = load_model()
    yres = species_result(4932, ypairs, [g.id for g in ymodel.genes])
    print(json.dumps(yres, indent=2))

    emodel = cobra.io.load_json_model(str(DATA / "iML1515.json"))
    egenes = json.load(open(RESULTS / "overrescue_audit_ecoli.json"))[
        "overrescued_genes"]
    epairs = pd.DataFrame(
        [{"gene": g, "backup": b}
         for g in egenes for b in sorted(isozyme_partners(emodel, g))])
    eres = species_result(511145, epairs, [g.id for g in emodel.genes])
    print(json.dumps(eres, indent=2))
    with open(RESULTS / "string_channels.json", "w") as fh:
        json.dump({"yeast": yres, "ecoli": eres}, fh, indent=2)


if __name__ == "__main__":
    main()
