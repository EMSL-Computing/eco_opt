from examples.iso_constraint import solve, REGRESSION_OBJECTIVE


def test_iso_constraint_satisfies_isoperimetric_constraint_and_regression():
    sol, sys = solve()
    assert sol.success

    # The isoperimetric constraint int_0^1 u^2 dt = 10 is enforced as a hard
    # equality on the auxiliary state x2(tf); verify IPOPT actually
    # satisfies it (not just reports success).
    assert abs(sol.x[-1, 1] - 10.0) < 1e-5

    # boundary conditions
    assert abs(sol.x[0, 0] - 1.0) < 1e-6
    assert abs(sol.x[0, 1] - 0.0) < 1e-6
    assert abs(sol.x[-1, 0] - 0.0) < 1e-5
    assert abs(sol.tf - 1.0) < 1e-6

    # No external reference value found for this exact problem instance;
    # regression-track eco_opt's own converged objective to catch silent
    # future regressions.
    assert abs(sol.objective - REGRESSION_OBJECTIVE) < 1e-4
