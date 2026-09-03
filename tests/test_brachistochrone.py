"""Analytic cycloid tests for the brachistochrone example."""
import math
import os
import sys as _sys

import numpy as np

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "brachistochrone"))

import brachistochrone as br


def test_cycloid_helpers_recover_half_cycloid():
    a, phi = br.cycloid_radius(br.X_F, br.Y_F)
    assert abs(a - 1.0) < 1e-12
    assert abs(phi - math.pi) < 1e-12
    assert abs(br.analytic_tf() - math.pi) < 1e-12
    assert abs(br.X_F - math.pi) < 1e-12
    assert abs(br.Y_F - 2.0) < 1e-12


def test_ocp_matches_analytic_cycloid():
    sol, sys = br.solve()
    assert sol.success, sol.return_status
    assert abs(sol.tf - br.ANALYTIC_TF) < 2e-3
    assert abs(sol.objective - br.ANALYTIC_TF) < 2e-3

    assert abs(sol.x[0, 0]) < 1e-8
    assert abs(sol.x[0, 1]) < 1e-8
    assert abs(sol.x[0, 2]) < 1e-8
    assert abs(sol.x[-1, 0] - br.X_F) < 1e-6
    assert abs(sol.x[-1, 1] - br.Y_F) < 1e-6

    xa, ya, va, tha = br.analytic_trajectories(sol.time)
    assert np.max(np.abs(sol.x[:, 0] - xa)) < 5e-3
    assert np.max(np.abs(sol.x[:, 1] - ya)) < 5e-3
    assert np.max(np.abs(sol.x[:, 2] - va)) < 2e-2

    # Energy: v^2 = 2 g y along a frictionless wire (skip t=0).
    g_y = 2.0 * br.G * sol.x[1:, 1]
    assert np.max(np.abs(sol.x[1:, 2] ** 2 - g_y)) < 5e-3

    # Heading tracks φ/2 after the v=0 instant.
    t_u = sol.time_intervals[:-1]
    tha_u = br.analytic_trajectories(t_u)[3]
    interior = t_u > 0.2
    assert np.mean(np.abs(sol.u[interior, 0] - tha_u[interior])) < 0.05
    assert sol.u[-1, 0] > 1.4  # near π/2
    assert np.all(sol.u[:, 0] >= -1e-6)
    assert np.all(sol.u[:, 0] <= math.pi + 1e-6)
