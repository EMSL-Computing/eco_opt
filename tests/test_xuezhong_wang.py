"""Regression + analytic-reference tests for the Xuezhong Wang examples."""
import os
import sys as _sys

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "xuezhong_wang"))

import bryson_denham_fixed_tf as ex1
import double_integrator_nonneg_control as ex2


def test_ex1_matches_analytic_lq_solution():
    sol, sys = ex1.solve(control_intervals=100, polynomial_degree=5)
    assert sol.success
    assert abs(sol.objective - ex1.ANALYTIC_OBJECTIVE) < 5e-3
    assert abs(sol.u[0, 0] - ex1.analytic_control(0.0)) < 5e-2
    assert abs(sol.u[-1, 0] - ex1.analytic_control(2.0)) < 5e-2
    assert abs(sol.x[0, 0]) < 1e-6
    assert abs(sol.x[0, 1]) < 1e-6
    assert abs(sol.x[-1, 0] - 5.0) < 1e-4
    assert abs(sol.x[-1, 1] - 2.0) < 1e-4


def test_ex2_trivial_zero_control_optimum():
    sol, sys = ex2.solve()
    assert sol.success
    # Constant-velocity coast (u==0) exactly satisfies all boundary
    # conditions, so this is the trivial global optimum, J*=0.
    assert abs(sol.objective) < 1e-6
    assert np.all(sol.u[:, 0] > -1e-4)  # control stays within a>=0
    assert abs(sol.x[0, 0] - 1.0) < 1e-6
    assert abs(sol.x[0, 1] - 1.0) < 1e-6
    assert abs(sol.x[-1, 0] - 3.0) < 1e-4
    assert sol.x[-1, 1] > -1e-4
