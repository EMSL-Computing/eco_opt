"""Minimum-energy double integrator, fixed final time (Xuezhong Wang, ex. 1).

Problem::

    minimize   int_0^2 0.5*a^2 dt
    subject to xdot = v, vdot = a
               x(0) = 0, v(0) = 0
               x(2) = 5, v(2) = 2

Analytic solution
------------------
This is a linear-quadratic (LQ) two-point boundary value problem: the
Euler-Lagrange / Pontryagin stationarity condition gives a costate that is
affine in t, hence an optimal control that is itself affine in time,
``a*(t) = alpha*t + beta``. Integrating twice::

    v(t) = alpha*t^2/2 + beta*t          (v(0)=0)
    x(t) = alpha*t^3/6 + beta*t^2/2      (x(0)=0)

and solving the two endpoint conditions ``v(2)=2``, ``x(2)=5`` for
``(alpha, beta)`` gives ``alpha=-4.5``, ``beta=5.5`` (see derivation in the
project session notes / re-derivable directly by solving the 2x2 linear
system above), so::

    a*(t) = -4.5*t + 5.5

and the analytic optimal objective is ``J* = 7.75`` (integrate
``0.5*(alpha*t+beta)^2`` from 0 to 2 with the above coefficients).
"""

import eco_opt
from eco_opt import Namespace


def trolley_model(t, x, u):
    position, speed = x[0], x[1]
    acceleration = u[0]
    y = Namespace(position=position, speed=speed, acceleration=acceleration)
    return [speed, acceleration], y


def build_ocp():
    sys = eco_opt.System(states=2, controls=1, model=trolley_model, name="trolley")
    time = sys.t
    trolley = sys.y

    ocp = eco_opt.Ocp()
    ocp.min([0.5 * trolley.acceleration ** 2])
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(time)),
            (0, "==", eco_opt.t_0(trolley.position)),
            (0, "==", eco_opt.t_0(trolley.speed)),
            (2, "==", eco_opt.t_f(time)),
            (5, "==", eco_opt.t_f(trolley.position)),
            (2, "==", eco_opt.t_f(trolley.speed)),
        ],
    )
    return ocp, sys


def solve(control_intervals=25, polynomial_degree=5):
    ocp, sys = build_ocp()
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        tf_guess=2.0,
    )
    return sol, sys


# Analytic reference values (see module docstring derivation).
ANALYTIC_OBJECTIVE = 7.75
ANALYTIC_ALPHA = -4.5
ANALYTIC_BETA = 5.5


def analytic_control(t):
    return ANALYTIC_ALPHA * t + ANALYTIC_BETA


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective:", sol.objective, "(analytic:", ANALYTIC_OBJECTIVE, ")")
    print("u(0):", sol.u[0, 0], "(analytic:", analytic_control(0.0), ")")
    print("u(tf):", sol.u[-1, 0], "(analytic:", analytic_control(2.0), ")")
