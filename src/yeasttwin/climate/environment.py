"""Combined-stress environments (pre-registration Section 5).

Axes and ranges are locked: T 28-42 C (1 C), ethanol 0-14% v/v (1%),
osmotic 0-1.5 M sorbitol-equivalent (0.1 M), nitrogen 5-100% of standard
(5%). Full grid 15 x 15 x 16 x 20 = 72,000. Combined-stress work uses a
5,000-environment Latin-hypercube sample (seed 20260926) with a stratified
70/30 DESIGN/TEST split (seed 424242, strata = temperature tercile x
ethanol tercile).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

GRID_SEED = 20260926
SPLIT_SEED = 424242
LHS_N = 5000
TEST_FRACTION = 0.30


@dataclass(frozen=True)
class Environment:
    temperature_c: float = 30.0
    ethanol_pct: float = 0.0
    osmotic_m: float = 0.0
    nitrogen_frac: float = 1.0


REFERENCE = Environment()


def grid_axes() -> dict[str, np.ndarray]:
    return {
        "temperature_c": np.arange(28.0, 42.0 + 1e-9, 1.0),
        "ethanol_pct": np.arange(0.0, 14.0 + 1e-9, 1.0),
        "osmotic_m": np.round(np.arange(0.0, 1.5 + 1e-9, 0.1), 10),
        "nitrogen_frac": np.round(np.arange(0.05, 1.0 + 1e-9, 0.05), 10),
    }


def full_grid() -> list[Environment]:
    ax = grid_axes()
    return [
        Environment(float(t), float(e), float(o), float(n))
        for t in ax["temperature_c"] for e in ax["ethanol_pct"]
        for o in ax["osmotic_m"] for n in ax["nitrogen_frac"]
    ]


def lhs_sample(n: int = LHS_N, seed: int = GRID_SEED) -> list[Environment]:
    """Latin-hypercube sample over the four locked ranges, reproducible by seed."""
    ax = grid_axes()
    rng = np.random.default_rng(seed)
    d = 4
    cut = np.linspace(0, 1, n + 1)
    u = rng.random((n, d))
    a, b = cut[:n], cut[1:n + 1]
    rd = a[:, None] + (b - a)[:, None] * u
    h = np.zeros_like(rd)
    for j in range(d):
        h[:, j] = rd[rng.permutation(n), j]
    keys = ["temperature_c", "ethanol_pct", "osmotic_m", "nitrogen_frac"]
    envs = []
    for i in range(n):
        vals = {}
        for j, k in enumerate(keys):
            lo, hi = float(ax[k][0]), float(ax[k][-1])
            vals[k] = round(lo + h[i, j] * (hi - lo), 6)
        envs.append(Environment(**vals))
    return envs


def design_test_split(envs: list[Environment] | None = None,
                      seed: int = SPLIT_SEED,
                      test_fraction: float = TEST_FRACTION
                      ) -> tuple[list[Environment], list[Environment]]:
    """Stratified 70/30 split (strata: T tercile x EtOH tercile), seeded."""
    envs = envs if envs is not None else lhs_sample()
    rng = np.random.default_rng(seed)
    ax = grid_axes()
    t_edges = np.quantile(ax["temperature_c"], [1 / 3, 2 / 3])
    e_edges = np.quantile(ax["ethanol_pct"], [1 / 3, 2 / 3])
    strata: dict[tuple[int, int], list[int]] = {}
    for i, env in enumerate(envs):
        s = (int(np.searchsorted(t_edges, env.temperature_c)),
             int(np.searchsorted(e_edges, env.ethanol_pct)))
        strata.setdefault(s, []).append(i)
    # largest-remainder apportionment: exact global 70/30 with strata honored
    target_test = int(round(len(envs) * test_fraction))
    raw = {s: len(idxs) * test_fraction for s, idxs in strata.items()}
    quotas = {s: int(v) for s, v in raw.items()}
    remainder = target_test - sum(quotas.values())
    for s in sorted(raw, key=lambda k: raw[k] - quotas[k], reverse=True)[:remainder]:
        quotas[s] += 1
    design_idx: list[int] = []
    test_idx: list[int] = []
    for s, idxs in strata.items():
        perm = rng.permutation(idxs)
        test_idx.extend(perm[:quotas[s]])
        design_idx.extend(perm[quotas[s]:])
    return ([envs[i] for i in sorted(design_idx)],
            [envs[i] for i in sorted(test_idx)])
