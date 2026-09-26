"""Climate-resilient fermentation: virtual S. cerevisiae cell under combined stress.

Redirect of MEGA27 item 5 (see REDIRECT.md). The scientific contract is
docs/PREREGISTRATION.md (SHA-256 b2235b022fc90d8791d524e9133507e688480651bd2c7551015f38662816dfc5).
Every parameter, range, seed, objective and gate in this package implements
that locked document; nothing here was tuned against outcome data.
"""
from .environment import Environment, REFERENCE, grid_axes, full_grid, lhs_sample, design_test_split
from .apply import applied
from .objectives import ethanol_yield, is_collapsed, COLLAPSE_FRACTION
