"""Minimum-time trolley stop, free terminal time (Lenhart & Workman,
ex. 20.4).

Problem::

    minimize   tf
    subject to xdot = v, vdot = a,  -1 <= a <= 1
               x(0) = x0 = 10, v(0) = 0
               x(tf) = 0, v(tf) >= 0

Analytic reference
-------------------
This is the classic minimum-time double-integrator "drive to rest at the
origin" problem. With |a|<=1 and v(0)=0, x(0)=x0>0, the time-optimal
control is bang-bang with exactly one switch: accelerate toward the
origin (a=-1) for the first half of the maneuver, then decelerate
(a=+1) for the second half, arriving at x=0 with v=0 exactly (an
interior optimum of the v(tf)>=0 constraint). By symmetry each half
covers a distance x0/2 starting/ending at rest under constant
|a|=1, so each half takes ``sqrt(x0)`` and the total minimum time is::

    tf* = 2*sqrt(x0)

For x0=10: ``tf* = 2*sqrt(10) ~= 6.32456``.
"""

import math

import casadi as cs

import eco_opt
from eco_opt import Namespace

X0 = 10.0


def trolley_model(t, x, u):
    position, speed = x[0], x[1]
    acceleration = u[0]
    y = Namespace(position=position, speed=speed, acceleration=acceleration)
    return [speed, acceleration], y


def build_ocp(x0=X0):
    sys = eco_opt.System(states=2, controls=1, model=trolley_model, name="trolley_mintime")
    time = sys.t
    trolley = sys.y

    ocp = eco_opt.Ocp()
    ocp.min([eco_opt.t_f(time)])
    ocp.st(
        system=sys,
        constraints=[
            (-1, "<=", trolley.acceleration, "<=", 1),
            (0, "==", eco_opt.t_0(time)),
            (x0, "==", eco_opt.t_0(trolley.position)),
            (0, "==", eco_opt.t_0(trolley.speed)),
            (0, "<=", eco_opt.t_f(time), "<=", cs.inf),
            (0, "==", eco_opt.t_f(trolley.position)),
            (0, "<=", eco_opt.t_f(trolley.speed), "<=", cs.inf),
        ],
    )
    return ocp, sys


def solve(control_intervals=50, polynomial_degree=3, x0=X0, tf_guess=None):
    ocp, sys = build_ocp(x0=x0)
    if tf_guess is None:
        tf_guess = analytic_tf(x0)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        tf_guess=tf_guess,
    )
    return sol, sys


def analytic_tf(x0=X0):
    return 2 * math.sqrt(x0)


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("tf:", sol.tf, "(analytic 2*sqrt(x0):", analytic_tf(), ")")
    print("u:", sol.u[:5, 0], "...", sol.u[-5:, 0])
