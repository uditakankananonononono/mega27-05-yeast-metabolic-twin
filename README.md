# mega27-05-yeast-metabolic-twin

A virtual *Saccharomyces cerevisiae* cell for metabolic engineering:
genome-scale FBA + graph deep learning on the consensus yeast model
(yeast-GEM: 4,105 reactions / 2,748 metabolites / 1,143 genes).

## Headline results
- **Published essentiality benchmark reproduced exactly** (acc 0.9015, MCC 0.532).
- **Locked-protocol CV**: no model beats the FBA rule on thresholded MCC
  (honest negative, preserved); the GCN **breaks the benchmark on ranking**:
  AUC 0.831 vs 0.697, paired delta +0.135, CI95 [0.114, 0.155].
- **Discovery: isozyme over-rescue** — FBA-missed essentials carry model
  isozyme backups 20.2% vs 1.5% (OR 16.2, p=3.8e-4); 19 genes named;
  reusable auditor: `python -m yeasttwin.overrescue`.
- **Discovery: Δsdh2Δach1** — obligatory succinate 0.0892 at 99.2% WT growth,
  +49% over the Raab-2010 quadruple with 2 deletions instead of 4.
  SDH2/SDH3 singles alone (0.067) also beat the quadruple.
- **Paper**: `paper/main.pdf` (50 pages, Times, 9 numbered equations with
  proofs, 5 figures, 8 real-data tables incl. full 1,107-gene tool output).

## Layout
- `src/yeasttwin/` — model loading, media, features, GCN/CNN learners,
  locked-gate CV harness, over-rescue auditor, strain-design scans
- `data/raw/` — SBML fixture, gold labels, sequences (all in-repo, hermetic)
- `results/` — every number in the paper, as CSV/JSON
- `tests/` — 76 hermetic test functions (no network)
- `scripts/` — scan runners, figure/table generators
- `paper/` — LaTeX source + compiled 50-page PDF

## Reproduce
```
pip install -r requirements.txt
PYTHONPATH=src python -m pytest tests/          # 76 test functions, hermetic
PYTHONPATH=src python -m yeasttwin.evaluate     # locked-gate CV
PYTHONPATH=src python -m yeasttwin.overrescue   # auditor
PYTHONPATH=src python scripts/run_strain_scan.py
```

## overrescue-audit (CLI tool)

Audits any genome-scale metabolic model for **isozyme over-rescue**:
essential genes the model calls viable only because an alternative-isozyme
GPR clause offers a backup that does not exist in vivo. Discovered in
yeast-GEM (OR 16.2), replicated in E. coli iML1515 (OR 11.5),
B. subtilis iYO844 (OR 22.7), and the older iMM904 yeast reconstruction.

```bash
./overrescue-audit --model model.json --essentials essentials.txt --out report.json
```

- `--model`: SBML (.xml) or JSON, COBRApy-readable
- `--essentials`: one essential gene ID per line (model gene IDs)
- report: Fisher odds ratio + p-value, missed/caught isozyme-backed
  fractions, and the named over-rescued genes

Verified: the CLI output on iML1515 + the Keio essential list
(`data/raw/keio_essential_baba2006_supTable6.xls`) reproduces
`results/overrescue_audit_ecoli.json` exactly.
