import sys
sys.path.insert(0, "src")
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update({"font.size": 10, "figure.dpi": 150})

# Fig 1: MCC + AUC distributions per model across locked folds
df = pd.read_csv("results/cv_fold_metrics.csv")
order = ["fba_rule", "logreg", "gcn", "cnn", "gcn_cnn"]
labels = ["FBA rule\n(published)", "Logistic\nregression", "GCN", "Seq CNN", "GCN+CNN"]
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
for ax, metric, title in zip(axes, ["mcc", "auc"],
                             ["MCC (locked gate)", "AUC-ROC (ranking)"]):
    data = [df[df.model == m][metric] for m in order]
    ax.boxplot(data, labels=labels, showmeans=True)
    ax.axhline(0.532 if metric == "mcc" else 0.697, ls="--", c="gray", lw=1)
    ax.set_title(title); ax.tick_params(axis="x", labelsize=8)
axes[0].set_ylabel("score across 15 folds")
fig.tight_layout(); fig.savefig("figures/fig_cv_metrics.png")

# Fig 2: succinate production envelope (recomputed from model)
from yeasttwin.model import load_model
from yeasttwin.strain_design import wt_production_envelope
env = wt_production_envelope(load_model(), points=10)
fig, ax = plt.subplots(figsize=(4.2, 3.4))
ax.plot(env["growth_frac"], env["succinate_flux"], "o-")
ax.set_xlabel("enforced growth (fraction of WT optimum)")
ax.set_ylabel("max succinate flux")
ax.set_title("Wild-type succinate production envelope")
fig.tight_layout(); fig.savefig("figures/fig_envelope.png")

# Fig 3: isozyme over-rescue
aud = json.load(open("results/overrescue_audit.json"))
fig, ax = plt.subplots(figsize=(4.2, 3.4))
ax.bar(["FBA-caught\nessentials (n=65)", "FBA-missed\nessentials (n=94)"],
       [aud["caught_isozyme_backed_frac"], aud["missed_isozyme_backed_frac"]],
       color=["#4477AA", "#CC6677"])
ax.set_ylabel("fraction with model isozyme backup")
ax.set_title("Isozyme over-rescue (OR=%.1f, p=%.1e)" % (aud["fisher_odds_ratio"], aud["fisher_p"]))
fig.tight_layout(); fig.savefig("figures/fig_overrescue.png")
print("figures done")
