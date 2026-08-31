import numpy as np

from examples.non_diff import solve, REGRESSION_OBJECTIVE


def test_non_diff_converges_and_regression():
    sol, sys = solve()
    # Nondifferentiable if_else forcing is a known-delicate case; the
    # requirement is that IPOPT actually reports Solve_Succeeded, not just
    # "ran without exception".
    assert sol.success
    assert sol.return_status == "Solve_Succeeded"

    assert abs(sol.t0 - 0.0) < 1e-9
    assert abs(sol.tf - 2.0) < 1e-6

    # x3(tf) (accumulated running cost) must be a small positive finite
    # number, not NaN/huge/negative.
    obj = sol.objective
    assert np.isfinite(obj)
    assert 0.0 < obj < 1000.0

    # transient response: x1/x2 should move away from zero due to the pulse
    assert sol.x[:, 0].max() - sol.x[:, 0].min() > 0.1
    assert sol.x[:, 1].max() - sol.x[:, 1].min() > 0.1

    # No .mat/external reference available; regression-track eco_opt's own
    # converged value (radau, degree=3, N=100, initial_guess={u: -5.0}).
    assert abs(obj - REGRESSION_OBJECTIVE) < 1e-2
