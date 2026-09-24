"""Train all models on the full labeled set, save per-gene scores, plot ROC."""
import sys
sys.path.insert(0, "src")
import numpy as np, pandas as pd, torch, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from yeasttwin.evaluate import (load_feature_table, load_graph, load_sequences)
from yeasttwin.labels import load_labels
from yeasttwin.ml import encode_sequences, train_cnn, train_gcn

labels = load_labels()
feat = load_feature_table()
genes = [g for g in feat.index if g in labels.index]
y = labels.loc[genes].to_numpy(bool)
x = feat.loc[genes].to_numpy(np.float32)
ggenes, a_hat = load_graph()
gpos = {g: i for i, g in enumerate(ggenes)}
order = np.array([gpos[g] for g in genes])
sc = StandardScaler().fit(x)
xg_full = np.zeros((len(ggenes), x.shape[1]), np.float32)
xg_full[order] = sc.transform(x)
y_full = np.zeros(len(ggenes), np.float32); y_full[order] = y.astype(np.float32)
tok = encode_sequences(load_sequences(ggenes))

s_fba = 1 - x[:, 0] * sc.scale_[0] - sc.mean_[0]
lr = LogisticRegression(max_iter=2000, class_weight="balanced").fit(sc.transform(x), y)
s_lr = lr.predict_proba(sc.transform(x))[:, 1]
s_gcn = train_gcn(torch.tensor(xg_full), a_hat, y_full, torch.tensor(order), epochs=300, seed=1)[order]
s_cnn = train_cnn(tok, y_full, order, epochs=25, seed=1)[order]

out = pd.DataFrame({"gene": genes, "essential": y.astype(int), "score_fba": s_fba,
                    "score_logreg": s_lr, "score_gcn": s_gcn, "score_cnn": s_cnn})
out = out.sort_values("score_gcn", ascending=False)
out.to_csv("results/full_model_gene_scores.csv", index=False)

fig, ax = plt.subplots(figsize=(4.6, 4.2))
for col, lab, c in [("score_fba", "FBA rule (AUC 0.697)", "#999999"),
                    ("score_logreg", "LogReg", "#4477AA"),
                    ("score_gcn", "GCN (AUC 0.831)", "#117733"),
                    ("score_cnn", "Seq CNN", "#CC6677")]:
    srt = out.sort_values(col, ascending=False)
    tps = np.cumsum(srt.essential) / y.sum()
    fps = np.cumsum(1 - srt.essential) / (~y).sum()
    ax.plot(np.r_[0, fps], np.r_[0, tps], label=lab, color=c)
ax.plot([0, 1], [0, 1], "k:", lw=0.8)
ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate")
ax.set_title("Essentiality ranking (full-data fit)")
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig("figures/fig_roc.png")
print("roc + scores done")
