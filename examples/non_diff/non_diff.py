"""Nondifferentiable system with a discontinuous pulse forcing term.

Exercises: model with no ``y`` signal tree; fixed final time (t_f(t)==2);
CasADi ``if_else``; radau collocation with polynomial degree 3; control
initial guess of -5.

No closed-form objective is available; this problem is numerically
delicate (nondifferentiable forcing). The converged objective is
regression-tracked.
"""

import casadi as cs

import eco_opt


def nondiff_model(t, x, u):
    x1, x2, x3 = x[0], x[1], x[2]
    d = 100.0 * (cs.if_else(t > 0.5, 1.0, 0.0) - cs.if_else(t > 0.6, 1.0, 0.0))
    dx1 = x2
    dx2 = -x1 - x2 + u[0] + d
    dx3 = 5 * x1 ** 2 + 2.5 * x2 ** 2 + 0.5 * u[0] ** 2
    return [dx1, dx2, dx3]  # no `y` signal tree


def build_ocp():
    sys = eco_opt.System(states=3, controls=1, model=nondiff_model, name="nondiff")
    time = sys.t
    x1, x2, x3 = sys.x[0], sys.x[1], sys.x[2]

    ocp = eco_opt.Ocp()
    ocp.min(eco_opt.t_f(x3))
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(time)),
            (0, "==", eco_opt.t_0(x1)),
            (0, "==", eco_opt.t_0(x2)),
            (0, "==", eco_opt.t_0(x3)),
            (2, "==", eco_opt.t_f(time)),
        ],
    )
    return ocp, sys


def solve(control_intervals=100, polynomial_degree=3):
    ocp, sys = build_ocp()
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="radau",
        initial_guess={sys.u[0]: -5.0},
        tf_guess=2.0,
    )
    return sol, sys


# Regression-tracked value: no external reference available; computed with
# eco_opt at control_intervals=100, polynomial_degree=3, radau collocation,
# initial_guess={u: -5.0}
REGRESSION_OBJECTIVE = 58.06949994233442


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective (x3(tf)):", sol.objective)
