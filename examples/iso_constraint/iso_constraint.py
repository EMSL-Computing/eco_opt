"""Isoperimetric constraint problem.

Exercises a model with no ``y`` signal tree, an auxiliary integral state
enforcing an isoperimetric constraint (int_0^1 u^2 dt = 10) as a fixed
box_final value, and a bare Lagrange objective term.

No external reference value found for this exact problem instance; we
report eco_opt's own converged objective and add it as a fixed-value
regression test to guard against silent regressions.
"""

import casadi as cs

import eco_opt


def iso_model(t, x, u):
    x1, x2 = x[0], x[1]
    return [-cs.sin(x1) + u[0], u[0] ** 2]  # no `y` signal tree


def build_ocp():
    sys = eco_opt.System(states=2, controls=1, model=iso_model, name="iso")
    time = sys.t
    x1, x2 = sys.x[0], sys.x[1]
    u = sys.u[0]

    ocp = eco_opt.Ocp()
    ocp.min([x1])  # LagrangeTerm(x1): integral of x1 dt over the horizon
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(time)),
            (1, "==", eco_opt.t_0(x1)),
            (0, "==", eco_opt.t_0(x2)),
            (-4, "<=", u, "<=", 4),
            (-10, "<=", x1, "<=", 10),
            (-float("inf"), "<=", x2, "<=", float("inf")),
            (1, "<=", eco_opt.t_f(time), "<=", 1),
            (0, "<=", eco_opt.t_f(x1), "<=", 0),
            (10, "==", eco_opt.t_f(x2)),
        ],
    )
    return ocp, sys


def solve(control_intervals=30, polynomial_degree=5):
    ocp, sys = build_ocp()
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        tf_guess=1.0,
    )
    return sol, sys


# Regression-tracked value: no external reference found for this exact
# instance; computed with eco_opt at control_intervals=30, degree=5,
# legendre collocation, default (1.0) initial guess.
REGRESSION_OBJECTIVE = -0.3751266963856673


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective (integral of x1):", sol.objective)
    print("x2(tf) (isoperimetric constraint, should be 10):", sol.x[-1, 1])
