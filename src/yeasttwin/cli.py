"""overrescue-audit: command-line isozyme over-rescue auditor.

Audits any genome-scale metabolic model (SBML or JSON, COBRApy-readable)
against a list of experimentally essential genes for the named failure
mode "isozyme over-rescue": essential genes the model calls viable
because an alternative-isozyme GPR clause offers a backup that does not
exist in vivo.

Usage:
    overrescue-audit --model model.json --essentials essentials.txt \
        [--growth-cutoff 1e-6] [--out report.json]

essentials.txt: one gene ID per line (must match model gene IDs).
Exit code 0 always; the JSON report carries the statistics.
"""
from __future__ import annotations

import argparse
import json
import sys

import cobra
import pandas as pd
from cobra.flux_analysis import single_gene_deletion

from .overrescue import overrescue_audit


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="overrescue-audit",
                                 description=__doc__.splitlines()[0])
    ap.add_argument("--model", required=True, help="SBML (.xml) or JSON model")
    ap.add_argument("--essentials", required=True,
                    help="text file: one essential gene ID per line")
    ap.add_argument("--growth-cutoff", type=float, default=1e-6,
                    help="KO/WT growth ratio below which a KO is inviable")
    ap.add_argument("--out", default=None, help="write JSON report here")
    args = ap.parse_args(argv)

    model = (cobra.io.load_json_model(args.model)
             if args.model.endswith(".json")
             else cobra.io.read_sbml_model(args.model))
    wt = model.slim_optimize()
    essentials = [l.strip() for l in open(args.essentials)
                  if l.strip() and not l.startswith("#")]
    in_model = [g for g in essentials if g in model.genes]
    ko = single_gene_deletion(model, gene_list=in_model, processes=1)
    ko = ko.dropna(subset=["growth"])
    ratios = pd.Series({next(iter(r["ids"])): float(r["growth"]) / wt
                        for _, r in ko.iterrows()})
    out = overrescue_audit(ratios, pd.Series(True, index=in_model),
                           model=model, viable_cutoff=args.growth_cutoff)
    out["model_file"] = args.model
    out["n_essentials_input"] = len(essentials)
    out["n_in_model"] = len(in_model)
    text = json.dumps(out, indent=2)
    print(text)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(text + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
