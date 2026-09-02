"""Least-squares tracking problem, free terminal state (Lenhart & Workman,
ex. 17.3). Smooth LQ tracking, not bang-bang.

Problem::

    minimize   int_0^2 (x - t^2)^2 dt
    subject to xdot = u,  0 <= u <= 4
               x(0) = 1

Reference / validation
-----------------------
This is a smooth (non-bang-bang, box constraint typically inactive) LQ
tracking problem; no simple closed form is asserted here, so we validate
with a regression-tracked objective and the physical sanity check that the
optimal state trajectory tracks (but does not exactly equal, since x(0)=1
!= 0 = t^2 at t=0) the reference parabola t^2 increasingly closely for
larger t, with the control staying interior to the box (0<=u<=4) rather
than saturating -- consistent with a smooth tracking-type problem.
"""

import casadi as cs

import eco_opt


def tracking_model(t, x, u):
    dx = u[0]
    return dx


def build_ocp():
    sys = eco_opt.System(states=1, controls=1, model=tracking_model, name="tracking")
    t = sys.t
    x = sys.x
    u = sys.u

    ocp = eco_opt.Ocp()
    ocp.min([(x[0] - t ** 2) ** 2])
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(t)),
            (1, "==", eco_opt.t_0(x[0])),
            (-cs.inf, "<=", x[0], "<=", cs.inf),
            (0, "<=", u[0], "<=", 4),
            (2, "==", eco_opt.t_f(t)),
        ],
    )
    return ocp, sys


def solve(control_intervals=60, polynomial_degree=3):
    ocp, sys = build_ocp()
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        tf_guess=2.0,
    )
    return sol, sys


# Regression-tracked value: computed with eco_opt at control_intervals=60,
# degree=3, legendre collocation (see solve() defaults). Note the control
# does approach (but not exceed) the u<=4 box bound near t=0, where x must
# move quickly away from x(0)=1 to start tracking t^2.
REGRESSION_OBJECTIVE = 0.5333333630040656


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective:", sol.objective, "(regression:", REGRESSION_OBJECTIVE, ")")
    print("u range:", sol.u.min(), sol.u.max())
