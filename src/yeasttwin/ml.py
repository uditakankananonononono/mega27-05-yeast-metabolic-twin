"""Learning layer: graph convolutional network over the metabolic gene graph
and a 1-D convolutional network over protein sequence, evaluated against the
FBA knockout rule under an identical gene-held-out cross-validation protocol.

GCN propagation (Kipf & Welling 2017):
    H^{(l+1)} = ReLU( D^{-1/2} (A + I) D^{-1/2} H^{(l)} W^{(l)} )
Training is transductive-masked: all node features are visible, but the loss
only uses labels of training-fold genes; test-fold labels are never read.
"""
from __future__ import annotations

import math

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as Fn

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i + 1 for i, a in enumerate(AA)}  # 0 = pad / unknown


# ----------------------------------------------------------------- metrics
def confusion(y: np.ndarray, p: np.ndarray) -> tuple[int, int, int, int]:
    y = np.asarray(y).astype(bool); p = np.asarray(p).astype(bool)
    return int((y & p).sum()), int((~y & p).sum()), int((y & ~p).sum()), int((~y & ~p).sum())


def mcc(y, p) -> float:
    tp, fp, fn, tn = confusion(y, p)
    d = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return (tp * tn - fp * fn) / math.sqrt(d) if d else 0.0


def auc_roc(y, s) -> float:
    """Mann-Whitney AUC with tie handling."""
    y = np.asarray(y).astype(bool); s = np.asarray(s, dtype=float)
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s)); ranks[order] = np.arange(1, len(s) + 1)
    # average ranks for ties
    _, inv, cnt = np.unique(s, return_inverse=True, return_counts=True)
    sums = np.bincount(inv, weights=ranks); ranks = sums[inv] / cnt[inv]
    npos, nneg = y.sum(), (~y).sum()
    if npos == 0 or nneg == 0:
        return float("nan")
    return float((ranks[y].sum() - npos * (npos + 1) / 2) / (npos * nneg))


def auc_pr(y, s) -> float:
    """Average precision."""
    y = np.asarray(y).astype(bool); s = np.asarray(s, dtype=float)
    order = np.argsort(-s, kind="mergesort"); y = y[order]
    tp = np.cumsum(y); prec = tp / np.arange(1, len(y) + 1)
    return float((prec * y).sum() / max(y.sum(), 1))


def best_threshold(y, s) -> float:
    """Threshold maximising MCC on (training) scores."""
    cand = np.unique(s)
    if len(cand) > 400:
        cand = np.quantile(s, np.linspace(0, 1, 401))
    best, bt = -2.0, 0.5
    for t in cand:
        m = mcc(y, s >= t)
        if m > best:
            best, bt = m, t
    return float(bt)


# ----------------------------------------------------------------- graph
def normalized_adjacency(adj_lists: list[list[int]], n: int) -> torch.Tensor:
    rows, cols = [], []
    for i, nb in enumerate(adj_lists):
        for j in nb:
            rows.append(i); cols.append(j)
    rows += list(range(n)); cols += list(range(n))
    idx = torch.tensor([rows, cols], dtype=torch.long)
    val = torch.ones(idx.shape[1])
    deg = torch.zeros(n).scatter_add_(0, idx[0], val)
    dinv = deg.pow(-0.5)
    val = dinv[idx[0]] * val * dinv[idx[1]]
    return torch.sparse_coo_tensor(idx, val, (n, n)).coalesce()


class GCN(nn.Module):
    def __init__(self, d_in: int, d_hid: int = 64, n_layers: int = 2, dropout: float = 0.3):
        super().__init__()
        dims = [d_in] + [d_hid] * n_layers
        self.layers = nn.ModuleList(nn.Linear(a, b) for a, b in zip(dims[:-1], dims[1:]))
        self.skip = nn.Linear(d_in, d_hid)
        self.out = nn.Linear(d_hid, 1)
        self.dropout = dropout

    def forward(self, x: torch.Tensor, a_hat: torch.Tensor) -> torch.Tensor:
        h = x
        for lin in self.layers:
            h = Fn.dropout(h, self.dropout, self.training)
            h = torch.relu(torch.sparse.mm(a_hat, lin(h)))
        h = h + torch.relu(self.skip(x))  # residual keeps node-intrinsic signal
        return self.out(h).squeeze(-1)


def train_gcn(x, a_hat, y, train_idx, epochs=300, lr=0.01, wd=5e-4, seed=0, **kw) -> np.ndarray:
    torch.manual_seed(seed)
    model = GCN(x.shape[1], **kw)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    yt = torch.tensor(y, dtype=torch.float32)
    tr = torch.tensor(train_idx, dtype=torch.long)
    pos = yt[tr].sum().clamp(min=1); neg = len(tr) - pos
    lossf = nn.BCEWithLogitsLoss(pos_weight=neg / pos)
    for _ in range(epochs):
        model.train(); opt.zero_grad()
        loss = lossf(model(x, a_hat)[tr], yt[tr])
        loss.backward(); opt.step()
    model.eval()
    with torch.no_grad():
        return torch.sigmoid(model(x, a_hat)).numpy()


# ----------------------------------------------------------------- sequence CNN
def encode_sequences(seqs: list[str], max_len: int = 800) -> torch.Tensor:
    out = torch.zeros(len(seqs), max_len, dtype=torch.long)
    for i, s in enumerate(seqs):
        ids = [AA_IDX.get(c, 0) for c in (s or "")[:max_len]]
        if ids:
            out[i, : len(ids)] = torch.tensor(ids)
    return out


class SeqCNN(nn.Module):
    """Embedding -> multi-width Conv1d -> global max pool -> linear."""

    def __init__(self, emb: int = 16, ch: int = 32, widths=(3, 7, 15), dropout: float = 0.3):
        super().__init__()
        self.emb = nn.Embedding(len(AA) + 1, emb, padding_idx=0)
        self.convs = nn.ModuleList(nn.Conv1d(emb, ch, w, padding=w // 2) for w in widths)
        self.drop = nn.Dropout(dropout)
        self.out = nn.Linear(ch * len(widths) + 1, 1)

    def forward(self, tok: torch.Tensor) -> torch.Tensor:
        e = self.emb(tok).transpose(1, 2)
        mask = (tok > 0).float().unsqueeze(1)
        feats = [(torch.relu(c(e)) * mask - (1 - mask) * 1e4).max(-1).values for c in self.convs]
        length = mask.sum(-1) / 800.0
        return self.out(self.drop(torch.cat(feats + [length], 1))).squeeze(-1)


def train_cnn(tok, y, train_idx, epochs=25, lr=3e-3, seed=0, batch=64) -> np.ndarray:
    torch.manual_seed(seed)
    model = SeqCNN()
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    yt = torch.tensor(y, dtype=torch.float32)
    tr = np.asarray(train_idx)
    pos = max(yt[tr].sum().item(), 1); neg = len(tr) - pos
    lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(neg / pos))
    rng = np.random.default_rng(seed)
    for _ in range(epochs):
        model.train(); rng.shuffle(tr)
        for b in range(0, len(tr), batch):
            ib = torch.tensor(tr[b:b + batch])
            opt.zero_grad(); lossf(model(tok[ib]), yt[ib]).backward(); opt.step()
    model.eval(); preds = []
    with torch.no_grad():
        for b in range(0, len(tok), 256):
            preds.append(torch.sigmoid(model(tok[b:b + 256])))
    return torch.cat(preds).numpy()
