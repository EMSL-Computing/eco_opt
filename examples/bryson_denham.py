"""Bryson-Denham problem: minimum-energy control with a state path constraint.

Direct translation of ``yoptimization/pages/examples/bryson/brysonDenham.md``
(Bryson & Ho, Applied Optimal Control, 1975). Exercises:
* a Lagrange (running-cost) objective term,
* a fully-fixed final time (t_f(time) == 1, i.e. lb=ub=1),
* a single-sided path constraint (position <= 1/9, no lower bound), which
  is the core `_unpack_tuple` fix for 3-tuples with '<='/'>=' only.

Reference / analytic solution
------------------------------
This is a classical benchmark with a known closed-form optimal cost. With
l = 1/9 the optimal trajectory "kisses" the boundary x=l over an interval
[t1, 1-t1] and the closed-form objective value is J* = 4 (see e.g. Bryson &
Ho 1975 sec 3.7, and numerous OC benchmark galleries -- PROPT's own solver
run on this exact problem reports f_k = 4.000021... for a comparable
discretization, see tomopt.com/docs/propt/tomlab_propt018.php). We use
J* = 4.0 as the reference value with a documented tolerance.
"""

import casadi as cs

import eco_opt
from eco_opt.system import Namespace


def trolley_model(t, x, u):
    position, speed = x[0], x[1]
    acceleration = u[0]
    dx = cs.vertcat(speed, acceleration)
    y = Namespace(position=position, speed=speed, acceleration=acceleration)
    return dx, y


def build_ocp():
    sys = eco_opt.System(states=2, controls=1, model=trolley_model, name="trolley")
    time = sys.t
    trolley = sys.y

    ocp = eco_opt.Ocp()
    ocp.min([0.5 * trolley.acceleration ** 2])
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(trolley.position)),
            (1, "==", eco_opt.t_0(trolley.speed)),
            (1, "==", eco_opt.t_f(time)),
            (0, "==", eco_opt.t_f(trolley.position)),
            (-1, "==", eco_opt.t_f(trolley.speed)),
            (1.0 / 9.0, ">=", trolley.position),
        ],
    )
    return ocp, sys


def solve(control_intervals=20, polynomial_degree=5):
    ocp, sys = build_ocp()
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        tf_guess=1.0,
    )
    return sol, sys


ANALYTIC_OBJECTIVE = 4.0
ANALYTIC_MAX_POSITION = 1.0 / 9.0


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective:", sol.objective, "(analytic:", ANALYTIC_OBJECTIVE, ")")
    print("max position:", sol.x[:, 0].max(), "(bound:", ANALYTIC_MAX_POSITION, ")")
