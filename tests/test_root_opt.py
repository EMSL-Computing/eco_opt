"""Regression + physical-sanity test for the root_opt plant C-N allocation
OCP (see examples/root_opt/README.md)."""
import os
import sys as _sys

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "root_opt"))

import plant_soil_nitrogen_allocation as m


def test_base_scenario_matches_regression_and_is_physical():
    sol, sys = m.solve()
    assert sol.success
    assert abs(sol.objective - m.REGRESSION_OBJECTIVE) < 1e-2

    # Boundary conditions
    for i, key in enumerate(["leafC", "rootC", "exuC", "SOMC", "micC", "rootN", "SOMN", "micN", "inorgN"]):
        assert abs(sol.x[0, i] - m.STATE_INIT[key]) < 1e-3

    # Mass-balance / physical sanity: all pools stay non-negative.
    assert np.all(sol.x >= -1e-6)
    # Leaf and root C should grow substantially over the 120-day season.
    assert sol.x[-1, 0] > sol.x[0, 0]
    assert sol.x[-1, 1] > sol.x[0, 1]
    # Mineral N should be net-depleted (plant uptake > mineralization) given
    # the fixed leaching_rate=0 base scenario.
    assert sol.x[-1, 8] < sol.x[0, 8]

    # Control bounds respected
    assert np.all(sol.u[:, 1] >= m.PARAMS["min_exudation"] - 1e-4)
    assert np.all(sol.u[:, 1] <= m.PARAMS["max_exudation"] + 1e-4)
    assert np.all(sol.u[:, 0] >= -1e-4)


def test_low_mineral_n_scenario_reduces_growth():
    """N-limitation sanity check: lowering initial mineral N should reduce
    cumulative leaf growth relative to the base scenario (physical
    consistency with the paper's N-limitation narrative)."""
    sol_base, _ = m.solve()
    state_low_n = dict(m.STATE_INIT)
    state_low_n["inorgN"] = 10.0
    sol_low_n, _ = m.solve(state_init=state_low_n)

    assert sol_low_n.success
    assert sol_low_n.objective < sol_base.objective
