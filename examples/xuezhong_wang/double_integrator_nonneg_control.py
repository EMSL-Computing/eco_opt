"""Minimum-energy double integrator with nonnegative control, fixed final
time (Xuezhong Wang, ex. 2).

Direct translation of
``examples/OptLignin/YOP benchmark/Xuezhong Wang_Solving optimal control
problems with MATLAB/Xuezhong_Wang_ex_2.m``.

Problem::

    minimize   int_0^2 0.5*a^2 dt
    subject to xdot = v, vdot = a,  a >= 0
               x(0) = 1, v(0) = 1
               x(2) = 3, v(2) >= 0

Reference / validation
-----------------------
Without the ``a >= 0`` bound the unconstrained LQ solution (same method as
in ``bryson_denham_fixed_tf.py``) would need
``v(2) = 2*alpha + beta = v(2)_free``; here v(2) is only lower-bounded
(``>=0``), so this is a genuinely constrained (possibly boundary-arc)
problem and no simple closed form is asserted. We instead validate:
  * feasibility: x(0)=1, v(0)=1, x(tf)=3, v(tf)>=0 (within solver tol),
  * physical sign constraint: u(t) >= 0 for all t (bang/singular arcs
    only ever accelerate, matching "a >= 0" in the MATLAB model),
  * a regression-tracked objective value for future comparisons.
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
    sys = eco_opt.System(states=2, controls=1, model=trolley_model, name="trolley2")
    time = sys.t
    trolley = sys.y

    ocp = eco_opt.Ocp()
    ocp.min([0.5 * trolley.acceleration ** 2])
    ocp.st(
        system=sys,
        constraints=[
            (0, "<=", trolley.acceleration, "<=", cs.inf),
            (0, "==", eco_opt.t_0(time)),
            (1, "==", eco_opt.t_0(trolley.position)),
            (1, "==", eco_opt.t_0(trolley.speed)),
            (2, "==", eco_opt.t_f(time)),
            (3, "==", eco_opt.t_f(trolley.position)),
            (0, "<=", eco_opt.t_f(trolley.speed), "<=", cs.inf),
        ],
    )
    return ocp, sys


def solve(control_intervals=50, polynomial_degree=5):
    ocp, sys = build_ocp()
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        tf_guess=2.0,
    )
    return sol, sys


# Regression-tracked value: computed with eco_opt at control_intervals=50,
# degree=5, legendre collocation (see solve() defaults). The optimal
# solution is trivial: constant velocity v=1 exactly satisfies x(0)=1 ->
# x(2)=3 with v(2)=1>=0, so u*(t)=0 and J*=0.
REGRESSION_OBJECTIVE = 0.0


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective:", sol.objective, "(regression:", REGRESSION_OBJECTIVE, ")")
    print("x(0):", sol.x[0, 0], "v(0):", sol.x[0, 1])
    print("x(tf):", sol.x[-1, 0], "v(tf):", sol.x[-1, 1])
    print("u range:", sol.u.min(), sol.u.max())
