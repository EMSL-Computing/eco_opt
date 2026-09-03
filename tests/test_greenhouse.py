"""Greenhouse OCP: tabulated light / outside temperature as interpolants."""
import os
import sys as _sys

import numpy as np

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "greenhouse"))

import greenhouse as gh
import greenhouse_from_file as gh_file


def test_greenhouse_regression_and_forcing():
    sol, sys = gh.solve()
    assert sol.success, sol.return_status
    assert abs(sol.t0 - 0.0) < 1e-12
    assert abs(sol.tf - gh.T_HORIZON) < 1e-10
    assert abs(sol.x[0, 0] - gh.X1_0) < 1e-12
    assert abs(sol.x[0, 1] - gh.X2_0) < 1e-10
    assert abs(sol.objective - gh.REGRESSION_OBJECTIVE) < 1e-6
    assert np.all(sol.u[:, 0] >= gh.U_MIN - 1e-8)
    assert np.all(sol.u[:, 0] <= gh.U_MAX + 1e-8)
    assert sol.u[:, 0].max() > 9.5
    assert sol.x[-1, 0] > 0.0

    # Interpolants recover the clipped sinusoids.
    t_mid = 12.0
    I_mid = float(gh.light(t_mid))
    T_mid = float(gh.outside_temp(t_mid))
    I_exact = max(
        0.0,
        800.0 * np.sin(4.0 * np.pi * t_mid / gh.T_HORIZON - 0.65 * np.pi),
    )
    T_exact = 15.0 + 10.0 * np.sin(4.0 * np.pi * t_mid / gh.T_HORIZON - 0.65 * np.pi)
    assert abs(I_mid - I_exact) < 1.0
    assert abs(T_mid - T_exact) < 0.05

    # Growth uses light: I is not identically zero, and is clipped at night.
    I_nodes = np.array([float(gh.light(ti)) for ti in sol.time])
    assert I_nodes.max() > 100.0
    assert I_nodes.min() < 1.0

    # Heating is used more while light is on (growth ∝ I * indoor T).
    I_u = np.array([float(gh.light(ti)) for ti in sol.time_intervals[:-1]])
    u = sol.u[:, 0]
    day = u[I_u > 100.0]
    night = u[I_u < 1.0]
    assert day.size and night.size
    assert day.mean() > night.mean()


def test_greenhouse_from_file_matches_inline():
    assert os.path.isfile(gh_file.FORCING_CSV)
    sol_a, _ = gh.solve()
    sol_b, _ = gh_file.solve()
    assert sol_b.success, sol_b.return_status
    assert abs(sol_a.objective - sol_b.objective) < 1e-8
    assert np.allclose(sol_a.x, sol_b.x, atol=1e-8)
    assert np.allclose(sol_a.u, sol_b.u, atol=1e-8)
    assert abs(float(gh.light(12.0)) - float(gh_file.light(12.0))) < 1e-6
    assert abs(float(gh.outside_temp(12.0)) - float(gh_file.outside_temp(12.0))) < 1e-6


def test_greenhouse_mesh_refinement_agrees():
    sol, _ = gh.solve(control_intervals=100, polynomial_degree=5)
    sol_c, _ = gh.solve(control_intervals=40, polynomial_degree=3)
    assert sol.success and sol_c.success
    rel = abs(sol_c.objective - sol.objective) / abs(sol.objective)
    assert rel < 0.05
