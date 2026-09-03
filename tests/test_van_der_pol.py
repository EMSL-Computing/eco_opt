"""Regression tests for the Van der Pol CasADi-collocation instance."""
import os
import sys as _sys

import numpy as np

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "van_der_pol"))

import van_der_pol as vdp


def test_vdp_matches_regression_and_bounds():
    sol, sys = vdp.solve()
    assert sol.success, sol.return_status
    assert abs(sol.tf - vdp.T_HORIZON) < 1e-10
    assert abs(sol.objective - vdp.REGRESSION_OBJECTIVE) < 1e-6
    assert abs(sol.x[0, 0] - vdp.X1_0) < 1e-10
    assert abs(sol.x[0, 1] - vdp.X2_0) < 1e-10
    assert np.all(sol.u[:, 0] >= vdp.U_MIN - 1e-8)
    assert np.all(sol.u[:, 0] <= vdp.U_MAX + 1e-8)
    # The quadratic cost drives the oscillator to the origin even without
    # a terminal-state constraint.
    assert abs(sol.x[-1, 0]) < 0.05
    assert abs(sol.x[-1, 1]) < 0.05
    # Upper bound is active (PMP saturation) on this instance.
    assert sol.u[:, 0].max() > 0.95


def test_vdp_mesh_refinement_agrees():
    """Coarse and default meshes should give the same J to a few percent."""
    sol_fine, _ = vdp.solve(control_intervals=40, polynomial_degree=4)
    sol_coarse, _ = vdp.solve(control_intervals=20, polynomial_degree=3)
    assert sol_fine.success and sol_coarse.success
    rel = abs(sol_coarse.objective - sol_fine.objective) / sol_fine.objective
    assert rel < 0.03
