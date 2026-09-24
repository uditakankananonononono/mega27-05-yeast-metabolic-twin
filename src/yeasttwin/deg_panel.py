"""DEG over-rescue panel: 3 new species + 4 cross-study replications.

Essentiality from DEG 15.2 (Database of Essential Genes, Luo et al 2021,
NAR 49:D677) bulk annotation table deg_annotation_p.csv - every row is one
essential gene call under an accession-level study (DEG1xxx, each with its
own PMID). Gene identifiers in DEG are original-study gene names or locus
tags; we map them to model genes through the organism's RefSeq GFF (NCBI
Datasets), which carries both current and old locus tags plus gene names.

New species: S. aureus N315/iSB619 (DEG1002, Ji 2001 Science),
S. Typhimurium LT2/STM_v1_0 (DEG1011, Knuth 2004 Mol Microbiol),
H. pylori 26695/iIT341 (DEG1008, Salama 2004 J Bacteriol).
Cross-study replications: E. coli/iML1515 vs Gerdes 2003 (DEG1018,
independent of Keio/Baba 2006), B. subtilis/iYO844 vs Kobayashi 2003
(DEG1001, independent of Koo 2017), M. tuberculosis/iNJ661 vs
Griffin 2011 Tn-seq (DEG1025) and Zhang 2012 (DEG1027), both independent
of the Sassetti 2003 calls used earlier.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import cobra
import pandas as pd
from cobra.flux_analysis import single_gene_deletion

from .overrescue import overrescue_audit

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "raw"
RESULTS = ROOT / "results"

# (model file, DEG accession, GFF file, organism, study, output suffix)
SPECS = [
    ("iSB619", "DEG1002", "gff_saureus_n315.gff", "Staphylococcus aureus N315",
     "Ji2001+Forsyth2002+Ko2006", "saureus_deg1002"),
    ("STM_v1_0", "DEG1011", "gff_salmonella_lt2.gff", "Salmonella typhimurium LT2",
     "Knuth2004", "styphimurium_deg1011"),
    ("iIT341", "DEG1008", "gff_hpylori_26695.gff", "Helicobacter pylori 26695",
     "Salama2004", "hpylori_deg1008"),
    ("iML1515", "DEG1018", "gff_ecoli_mg1655.gff", "Escherichia coli MG1655",
     "Gerdes2003", "ecoli_gerdes_deg1018"),
    ("iYO844", "DEG1001", "gff_bsubtilis_168.gff", "Bacillus subtilis 168",
     "Kobayashi2003", "bsubtilis_kobayashi_deg1001"),
    ("iNJ661", "DEG1025", "gff_mtb_h37rv.gff",
     "Mycobacterium tuberculosis H37Rv", "Griffin2011Tnseq", "mtb_griffin_deg1025"),
    ("iNJ661", "DEG1027", "gff_mtb_h37rv.gff",
     "Mycobacterium tuberculosis H37Rv", "Zhang2012", "mtb_zhang_deg1027"),
]


def parse_gff(path: Path):
    """name-lower -> locus tag, plus every locus-tag variant (current+old)."""
    name2tag, tags = {}, set()
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("gene", "CDS"):
            continue
        attrs = dict(re.findall(r"([^=;]+)=([^;]*)", f[8]))
        name = attrs.get("gene") or attrs.get("Name")
        variants = [attrs.get("locus_tag")]
        if attrs.get("old_locus_tag"):
            variants += attrs["old_locus_tag"].replace("%2C", ",").split(",")
        good = [t for t in variants if t]
        for tag in good:
            tags.add(tag)
        if name:
            name2tag.setdefault(name.lower(), set()).update(good)
    return name2tag, tags


def deg_essentials(acc: str) -> list[str]:
    out = []
    with open(DATA / "deg_annotation_p.csv") as fh:
        for r in csv.reader(fh, delimiter=";", quotechar='"'):
            if len(r) > 8 and r[0] == acc:
                out.append(r[2].strip())
    return out


def map_to_loci(genes: list[str], name2tag, tags):
    """Each essential gene maps to ALL its locus-tag variants, so the model
    intersection matches whichever tag generation the model uses."""
    mapped, unmapped = [], []
    for g in genes:
        if g in tags:
            mapped.append(g)
            continue
        hits = set()
        for part in re.split(r"[/;]", g):
            p = part.strip().lower()
            if p in name2tag:
                hits.update(name2tag[p])
        if hits:
            mapped.extend(sorted(hits))
        else:
            unmapped.append(g)
    return mapped, unmapped


def run_spec(spec) -> dict:
    model_id, acc, gff, organism, study, suffix = spec
    model = cobra.io.load_json_model(str(DATA / f"{model_id}.json"))
    name2tag, tags = parse_gff(DATA / gff)
    deg = deg_essentials(acc)
    mapped, unmapped = map_to_loci(deg, name2tag, tags)
    model_genes = {g.id for g in model.genes}
    in_model = sorted({g for g in mapped if g in model_genes})
    wt = model.slim_optimize()
    if not wt or wt <= 1e-6:
        return {"species": f"{organism} ({model_id}) vs {study} ({acc})",
                "error": "model biomass blocked under default medium",
                "wt": float(wt or 0.0)}
    ko = single_gene_deletion(model, gene_list=in_model, processes=1)
    ko = ko.dropna(subset=["growth"])
    ratios = pd.Series({next(iter(r["ids"])): float(r["growth"]) / wt
                        for _, r in ko.iterrows()})
    out = overrescue_audit(ratios, pd.Series(True, index=in_model), model=model)
    out["species"] = f"{organism} ({model_id}) vs {study} ({acc})"
    out["deg_accession"] = acc
    out["study"] = study
    out["n_deg_essentials"] = len(deg)
    out["n_mapped"] = len(mapped)
    out["n_unmapped"] = len(unmapped)
    out["n_in_model"] = len(in_model)
    with open(RESULTS / f"overrescue_audit_{suffix}.json", "w") as fh:
        json.dump(out, fh, indent=2)
    return out


def main():
    rows = []
    for spec in SPECS:
        out = run_spec(spec)
        print(f"{spec[5]}: OR {out.get('fisher_odds_ratio', float('nan')):.1f} "
              f"p={out.get('fisher_p', float('nan')):.2e} "
              f"overrescued={len(out.get('overrescued_genes', []))} "
              f"mapped={out.get('n_mapped')}/{out.get('n_deg_essentials')} "
              f"in_model={out.get('n_in_model')}")
        rows.append({k: (v if not isinstance(v, list) else len(v))
                     for k, v in out.items()})
    pd.DataFrame(rows).to_csv(RESULTS / "deg_panel_summary.csv", index=False)


if __name__ == "__main__":
    main()
