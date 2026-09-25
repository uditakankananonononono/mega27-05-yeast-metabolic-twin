"""statsmodels multiple-testing lane: FDR re-analysis of the over-rescue panel.

The 12 species/model over-rescue audits each report a Fisher exact p-value
for isozyme-backup enrichment in FBA-missed essentials. This lane applies
Benjamini-Hochberg FDR control (statsmodels multipletests) across the whole
panel so claims survive panel-wide multiple testing, and writes the adjusted
table committed under results/.
"""
import glob
import json
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent.parent.parent


def main():
    rows = []
    for p in sorted(glob.glob(str(ROOT / "results" / "overrescue_audit*.json"))):
        d = json.loads(Path(p).read_text())
        name = Path(p).stem.replace("overrescue_audit", "").strip("_") or "yeast"
        rows.append({"audit": name,
                     "n_missed": d.get("n_missed"),
                     "missed_backed_frac": d.get("missed_isozyme_backed_frac"),
                     "caught_backed_frac": d.get("caught_isozyme_backed_frac"),
                     "fisher_odds_ratio": d.get("fisher_odds_ratio"),
                     "fisher_p": d.get("fisher_p")})
    df = pd.DataFrame(rows).dropna(subset=["fisher_p"])
    rej, padj, _, _ = multipletests(df.fisher_p, alpha=0.05, method="fdr_bh")
    df["padj_bh"] = padj
    df["significant_fdr05"] = rej
    df.to_csv(ROOT / "results" / "statsmodels_fdr.csv", index=False)
    summary = {
        "panel_size": int(len(df)),
        "n_significant_raw": int((df.fisher_p < 0.05).sum()),
        "n_significant_fdr05": int(rej.sum()),
        "yeast_discovery_padj": float(df.loc[df.audit == "yeast", "padj_bh"].iloc[0]),
        "method": "Benjamini-Hochberg via statsmodels.stats.multitest.multipletests",
    }
    (ROOT / "results" / "statsmodels_fdr_summary.json").write_text(
        json.dumps(summary, indent=2))
    print(summary)


if __name__ == "__main__":
    main()
