import os
import sys as _sys

import numpy as np

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "bang_free_tf"))

from bang_free_tf import solve, ANALYTIC_TF, ANALYTIC_SWITCH_TIME


def test_bang_free_tf_matches_analytic_min_time():
    sol, sys = solve()
    assert sol.success

    # Analytic bang-bang double-integrator minimum time (see module
    # docstring derivation): tf = 30.0 exactly.
    rel_err = abs(sol.tf - ANALYTIC_TF) / ANALYTIC_TF
    assert rel_err < 1e-2, f"tf={sol.tf} vs analytic {ANALYTIC_TF}, rel_err={rel_err}"

    # boundary conditions
    assert abs(sol.x[0, 0] - 0.0) < 1e-6
    assert abs(sol.x[0, 1] - 0.0) < 1e-6
    assert abs(sol.x[-1, 0] - 300.0) < 1e-3
    assert abs(sol.x[-1, 1] - 0.0) < 1e-3

    # bang-bang control structure: u ~ +1 before the switch, u ~ -2 after
    u = sol.u[:, 0]
    t_mid = sol.time_intervals[:-1]
    before = u[t_mid < ANALYTIC_SWITCH_TIME - 2.0]
    after = u[t_mid > ANALYTIC_SWITCH_TIME + 2.0]
    assert before.size and after.size
    assert np.allclose(before, 1.0, atol=0.05)
    assert np.allclose(after, -2.0, atol=0.05)

    # single switch located near the analytic switch time
    switch_idx = np.argmax(np.abs(np.diff(u)) > 1.0)
    switch_time = t_mid[switch_idx]
    assert abs(switch_time - ANALYTIC_SWITCH_TIME) < 2.0
