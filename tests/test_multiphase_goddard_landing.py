"""Two-phase Goddard landing (Yop goddardLanding example)."""
import os
import sys as _sys

import numpy as np

HERE = os.path.dirname(__file__)
_sys.path.insert(
    0, os.path.join(HERE, "..", "examples", "multiphase_goddard_landing")
)

import goddard_landing as gl


def test_goddard_landing_soft_touchdown_and_regression():
    sol, systems = gl.solve()
    assert sol.success, sol.return_status
    assert len(sol.phases) == 2
    assert abs(sol.t0) < 1e-8
    assert sol.switch_times[0] > 1.0
    assert sol.tf > sol.switch_times[0] + 1.0

    # Soft landing.
    assert abs(sol.x[-1, 0]) < 1e-3  # v(T)
    assert abs(sol.x[-1, 1]) < 1e-3  # h(T)
    assert sol.x[-1, 2] >= gl.MF - 1e-3

    # Bounds.
    assert np.min(sol.x[:, 1]) >= -1e-4
    assert np.min(sol.x[:, 2]) >= gl.MF - 1e-3
    assert np.max(sol.x[:, 2]) <= gl.M0 + 1e-6
    assert np.min(sol.u[:, 0]) >= -1e-6
    assert np.max(sol.u[:, 0]) <= gl.FM + 1e-5

    # Junction continuity.
    assert abs(sol.phases[0].tf - sol.phases[1].t0) < 1e-6
    assert np.allclose(sol.phases[0].x[-1], sol.phases[1].x[0], atol=1e-4)

    # Peak height is the phase-1 Mayer term; must be well above the pad.
    assert sol.objective > 1e3
    assert abs(sol.objective / gl.REGRESSION_PEAK_HEIGHT - 1.0) < 0.05
    assert abs(sol.tf / gl.REGRESSION_TF - 1.0) < 0.05


def test_goddard_landing_optimize_matches_solve():
    sol_s, _ = gl.solve()
    sol_o, _ = gl.optimize()
    assert sol_s.success and sol_o.success
    assert abs(sol_s.objective - sol_o.objective) / abs(sol_s.objective) < 1e-6
    assert abs(sol_s.tf - sol_o.tf) / abs(sol_s.tf) < 1e-6
    assert abs(sol_o.x[-1, 0]) < 1e-3
    assert abs(sol_o.x[-1, 1]) < 1e-3
