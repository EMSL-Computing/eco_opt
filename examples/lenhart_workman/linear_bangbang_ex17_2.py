"""Linear objective, bang-bang control (Lenhart & Workman, ex. 17.2).

Direct translation of
``examples/OptLignin/YOP benchmark/Examples from Lenhart and
Workman_Optimal control applied to biological models/
BangBangControl_lenhart_ex17_2.m``.

Problem::

    maximize   int_0^2 (2*x - 3*u) dt
    subject to xdot = x + u,  0 <= u <= 2
               x(0) = 5

Analytic (Pontryagin) reference
--------------------------------
Hamiltonian H = 2x - 3u + lambda*(x+u) = 2x + lambda*x + u*(lambda-3).
Costate: lambda' = -dH/dx = -(2+lambda), with free x(tf) (no terminal
state penalty) giving the transversality condition lambda(tf) = 0. Solving
lambda' = -(2+lambda), lambda(2) = 0 gives::

    lambda(t) = -2 + 2*exp(2-t)

The switching function is lambda(t) - 3; since u enters H affinely,
Pontryagin's maximum principle gives bang-bang control::

    u*(t) = 2   if lambda(t) - 3 > 0
    u*(t) = 0   if lambda(t) - 3 < 0

Evaluating: lambda(t) - 3 = -5 + 2*exp(2-t), which is positive for
t < 2 - ln(5/2) approx 1.083 and negative afterwards, i.e. a single
switch from u=2 (t<t_switch) to u=0 (t>t_switch) at
``t_switch = 2 - ln(2.5) ~= 1.0837``.
"""

import math

import casadi as cs

import eco_opt


def linear_bangbang_model(t, x, u):
    dx = x + u
    return dx


def build_ocp():
    sys = eco_opt.System(states=1, controls=1, model=linear_bangbang_model, name="bangbang")
    x = sys.x
    u = sys.u

    ocp = eco_opt.Ocp()
    ocp.max([2 * x[0] - 3 * u[0]])
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(sys.t)),
            (5, "==", eco_opt.t_0(x[0])),
            (-cs.inf, "<=", x[0], "<=", cs.inf),
            (0, "<=", u[0], "<=", 2),
            (2, "==", eco_opt.t_f(sys.t)),
        ],
    )
    return ocp, sys


def solve(control_intervals=30, polynomial_degree=3):
    ocp, sys = build_ocp()
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        tf_guess=2.0,
    )
    return sol, sys


# Analytic switching time (see module docstring derivation).
ANALYTIC_SWITCH_TIME = 2 - math.log(2.5)


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("analytic switch time:", ANALYTIC_SWITCH_TIME)
    print("u:", sol.u[:, 0])
