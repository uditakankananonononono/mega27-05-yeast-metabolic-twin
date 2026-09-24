import numpy as np

from yeasttwin.evaluate import make_folds


def test_make_folds_stratified_and_deterministic():
    y = np.array([0]*90 + [1]*10)
    f1 = make_folds(y, 5, 2, 123)
    f2 = make_folds(y, 5, 2, 123)
    assert (f1 == f2).all()
    for r in range(2):
        for k in range(5):
            assert (y[f1[:, r] == k].sum() == 2)
