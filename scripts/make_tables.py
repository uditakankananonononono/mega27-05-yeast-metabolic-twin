"""Generate LaTeX tables from real result files (reproducible)."""
import json
import pandas as pd

df = pd.read_csv("results/cv_fold_metrics.csv")
piv = df.pivot_table(index=["repeat", "fold"], columns="model", values=["mcc", "auc"])
with open("paper/tab_cv_folds.tex", "w") as f:
    f.write("\\begin{longtable}{cc rrrrr rrrrr}\n")
    f.write("\\caption{Per-fold locked-protocol results.}\\label{tab:folds}\\\\\n")
    f.write("\\toprule\n & & \\multicolumn{5}{c}{MCC} & \\multicolumn{5}{c}{AUC-ROC}\\\\\n")
    f.write("repeat & fold & FBA & LR & GCN & CNN & ens & FBA & LR & GCN & CNN & ens\\\\\n\\midrule\n")
    for (r, k), row in piv.iterrows():
        cells = [f"{row[(m, mod)]:.3f}" for m in ["mcc", "auc"]
                 for mod in ["fba_rule", "logreg", "gcn", "cnn", "gcn_cnn"]]
        f.write(f"{r} & {k} & " + " & ".join(cells) + "\\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n")

aud = json.load(open("results/overrescue_audit.json"))
genes = pd.read_csv("results/overrescued_genes_annotated.csv")
with open("paper/tab_overrescued.tex", "w") as f:
    f.write("\\begin{longtable}{ll}\n\\caption{The 19 isozyme over-rescued essential genes.}\\label{tab:overrescued}\\\\\n\\toprule\nORF & annotation\\\\\n\\midrule\n")
    for _, r in genes.iterrows():
        ann = r["annotation"].replace("&", "\\&").replace("%", "\\%")
        f.write(f"{r['gene']} & {ann}\\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n")
print("tables written")
