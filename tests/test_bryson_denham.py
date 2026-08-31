from examples.bryson_denham import solve, ANALYTIC_OBJECTIVE, ANALYTIC_MAX_POSITION


def test_bryson_denham_matches_analytic_solution():
    sol, sys = solve()
    assert sol.success
    assert abs(sol.tf - 1.0) < 1e-6  # fully fixed final time

    rel_err = abs(sol.objective - ANALYTIC_OBJECTIVE) / ANALYTIC_OBJECTIVE
    assert rel_err < 5e-3, f"objective={sol.objective} vs analytic {ANALYTIC_OBJECTIVE}"

    # single-sided path constraint: position never exceeds 1/9 (to solver
    # tolerance) and gets close to the bound (the "kissing arc").
    max_pos = sol.x[:, 0].max()
    assert max_pos <= ANALYTIC_MAX_POSITION + 1e-3
    assert abs(max_pos - ANALYTIC_MAX_POSITION) < 1e-3

    # boundary conditions
    assert abs(sol.x[0, 0] - 0.0) < 1e-6
    assert abs(sol.x[0, 1] - 1.0) < 1e-6
    assert abs(sol.x[-1, 0] - 0.0) < 1e-4
    assert abs(sol.x[-1, 1] - (-1.0)) < 1e-4
