import numpy as np

from yeasttwin.ml import auc_pr, auc_roc, best_threshold, confusion, mcc


def test_confusion_and_mcc_perfect():
    y = np.array([1, 1, 0, 0, 1, 0])
    p = np.array([1, 1, 0, 0, 0, 0])
    assert confusion(y, p) == (2, 0, 1, 3)
    assert mcc(y, y) == 1.0
    assert mcc(y, 1 - y) == -1.0


def test_mcc_published_benchmark_reproduction():
    y = np.array([1]*159 + [0]*948)
    p = np.array([1]*65 + [0]*94 + [1]*15 + [0]*933)
    assert abs(mcc(y, p) - 0.532) < 5e-4


def test_auc_roc_edge_cases():
    y = np.array([0, 0, 1, 1])
    s = np.array([0.1, 0.4, 0.35, 0.8])
    assert abs(auc_roc(y, s) - 0.75) < 1e-9
    assert auc_roc(y, y) == 1.0


def test_auc_pr_perfect():
    y = np.array([0, 1, 0, 1])
    s = np.array([0.1, 0.9, 0.2, 0.8])
    assert auc_pr(y, s) == 1.0


def test_best_threshold_recovers_split():
    rng = np.random.default_rng(0)
    y = rng.integers(0, 2, 200)
    s = y * 0.6 + rng.normal(0, 0.2, 200)
    t = best_threshold(y, s)
    assert mcc(y, s >= t) > 0.5
