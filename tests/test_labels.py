from yeasttwin.labels import load_labels


def test_label_convention_matches_sgd_and_benchmark():
    lab = load_labels()
    assert lab.sum() == 159
    assert len(lab) == 1107
