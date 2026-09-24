"""KEGG pathway enrichment of bacterial over-rescued sets (hypergeometric).

g:Profiler's organism set does not cover the bacterial species (checked
against its 1,113-organism list: only ecoligprjeb31222), so bacterial
enrichment runs directly on KEGG REST: per species, model genes are
mapped to KEGG pathways (link/pathway/<org>), and each pathway is tested
for over-representation of the over-rescued set among pathway-annotated
model genes with the hypergeometric survival function

    P(X >= a) = sum_{k=a}^{n} C(K,k) C(N-K,n-k) / C(N,n)      [eq. 12]

Bonferroni across tested pathways. Background = model genes with any
KEGG pathway annotation (metabolic scope of the audit).
"""
from __future__ import annotations

import json
import re
import time
import urllib.request
from pathlib import Path

import cobra
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw"
RESULTS = ROOT / "results"

# species: (kegg org, model json, audit json, gff for RS->old tag or None)
SPECIES = [
    ("ecoli", "eco", "iML1515", "overrescue_audit_ecoli.json", None),
    ("bsubtilis", "bsu", "iYO844", "overrescue_audit_bsubtilis.json", None),
    ("mtb", "mtu", "iNJ661", "overrescue_audit_mtb.json", None),
    ("saureus", "sau", "iSB619", "overrescue_audit_saureus_deg1002.json",
     "gff_saureus_n315.gff"),
    ("styphimurium", "stm", "STM_v1_0",
     "overrescue_audit_styphimurium_deg1011.json", None),
]


def kegg_get(path: str) -> str:
    with urllib.request.urlopen("https://rest.kegg.jp/" + path,
                                timeout=60) as r:
        return r.read().decode()


def gff_rs_to_old(path: Path) -> dict:
    m = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        attrs = dict(re.findall(r"([^=;]+)=([^;]*)", f[8]))
        cur = attrs.get("locus_tag")
        olds = (attrs.get("old_locus_tag") or "").replace("%2C", ",").split(",")
        for o in olds:
            if o:
                m.setdefault(cur, o)
    return m


def main():
    summary = {}
    for tag, org, model_id, audit, gff in SPECIES:
        model = cobra.io.load_json_model(str(DATA / f"{model_id}.json"))
        over = json.load(open(RESULTS / audit))["overrescued_genes"]
        if gff:  # S. aureus model uses SA_RS; KEGG uses SA old tags
            rs2old = gff_rs_to_old(DATA / gff)
            over_k = [rs2old.get(g, g) for g in over]
            model_k = [rs2old.get(g.id, g.id) for g in model.genes]
        else:
            over_k, model_k = over, [g.id for g in model.genes]
        link = kegg_get(f"link/pathway/{org}")
        gene2pws = {}
        for ln in link.strip().split("\n"):
            g, pw = ln.split("\t")
            gene2pws.setdefault(g.split(":", 1)[1], set()).add(
                pw.split(":", 1)[1])
        bg = {g for g in model_k if g in gene2pws}
        oset = {g for g in over_k if g in bg}
        pws = sorted({p for g in bg for p in gene2pws[g]})
        N, n = len(bg), len(oset)
        rows = []
        for pw in pws:
            members = {g for g in bg if pw in gene2pws[g]}
            K, a = len(members), len(members & oset)
            if a == 0:
                continue
            p = float(hypergeom.sf(a - 1, N, K, n))
            rows.append({"pathway": pw, "overlap": a, "size": K, "p": p})
        df = pd.DataFrame(rows)
        if len(df):
            df["p_bonf"] = (df.p * len(pws)).clip(upper=1.0)
            df = df.sort_values("p")
            df.to_csv(RESULTS / f"kegg_enrich_{tag}.csv", index=False)
        sig = df[df.p_bonf < 0.05] if len(df) else df
        summary[tag] = {"kegg_org": org, "n_overrescued": len(over),
                        "n_background": N, "n_over_in_bg": n,
                        "n_pathways_tested": len(pws),
                        "n_significant": int(len(sig)),
                        "top": df.head(5).to_dict("records") if len(df) else []}
        print(tag, "bg", N, "over", n, "sig", len(sig),
              sig.head(3).pathway.tolist() if len(sig) else "-")
        time.sleep(0.3)
    with open(RESULTS / "kegg_enrich_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
