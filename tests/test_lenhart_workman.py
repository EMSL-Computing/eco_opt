"""Regression + analytic-reference tests for the Lenhart & Workman
examples."""
import os
import sys as _sys

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "lenhart_workman"))

import linear_bangbang_ex17_2 as ex172
import lq_tracking_ex17_3 as ex173
import min_time_trolley_ex20_4 as ex204


def test_ex17_2_bangbang_switch_matches_analytic():
    sol, sys = ex172.solve()
    assert sol.success
    t = sol.time_intervals[:-1]
    u = sol.u[:, 0]
    switch_idx = np.argmax(u < 1.0)  # first interval where u has dropped
    switch_time = t[switch_idx]
    assert abs(switch_time - ex172.ANALYTIC_SWITCH_TIME) < 0.1
    assert u[0] > 1.9  # u=2 (max) before the switch
    assert u[-1] < 0.1  # u=0 after the switch


def test_ex17_3_tracking_regression_and_box_bounds():
    sol, sys = ex173.solve()
    assert sol.success
    assert abs(sol.objective - ex173.REGRESSION_OBJECTIVE) < 1e-3
    assert np.all(sol.u[:, 0] >= -1e-4)
    assert np.all(sol.u[:, 0] <= 4 + 1e-4)
    assert abs(sol.x[0, 0] - 1.0) < 1e-6


def test_ex20_4_min_time_matches_analytic():
    sol, sys = ex204.solve()
    assert sol.success
    assert abs(sol.tf - ex204.analytic_tf()) < 1e-2
    assert abs(sol.x[0, 0] - 10.0) < 1e-4
    assert abs(sol.x[-1, 0]) < 1e-3
    assert sol.x[-1, 1] > -1e-3
    # bang-bang: control near -1 for early time, near +1 for late time
    assert sol.u[0, 0] < -0.9
    assert sol.u[-1, 0] > 0.9
