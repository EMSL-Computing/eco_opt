"""Bang-bang minimum-time double integrator, free final time.

The model returns a bare ``dx`` (no Namespace).

Analytic reference (asymmetric bang-bang double integrator):
    xdot1 = x2, xdot2 = u,  u in [-2, 1]
    x(0) = (0, 0), x(tf) = (300, 0), minimize tf.

For a double integrator driven from rest to rest over distance D with an
accelerate phase at u=u1 followed by a decelerate phase at u=-u2 (u1, u2 > 0),
the standard bang-bang derivation gives:
    t1 (switch time): u1 * t1 = v1 (velocity at switch)
    t2 = v1 / u2  (decel phase duration, ends at v=0)
    distance = u1*t1^2/2 + v1^2/(2*u2) = (t1^2/2) * u1 * (1 + u1/u2)
    T = t1 + t2 = t1 * (1 + u1/u2)

With u1=1, u2=2, D=300:
    (t1^2/2)*1*(1.5) = 300  =>  t1^2 = 400  =>  t1 = 20
    T = t1 * 1.5 = 30

So the analytic minimum time is tf = 30, with a single switch from u=+1 to
u=-2 at t=20.

Run directly:  python examples/bang_free_tf/bang_free_tf.py
"""

import casadi as cs

import eco_opt


def bang_model(t, x, u):
    dx = cs.vertcat(x[1], u[0])
    return dx  # no `y` signal tree -- exercises System's optional-y support


def build_ocp():
    sys = eco_opt.System(states=2, controls=1, model=bang_model, name="bang")
    time = sys.t
    x1, x2 = sys.x[0], sys.x[1]

    ocp = eco_opt.Ocp()
    ocp.min(eco_opt.t_f(time))
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(time)),
            (0, "==", eco_opt.t_0(x1)),
            (0, "==", eco_opt.t_0(x2)),
            (-float("inf"), "<=", x1, "<=", float("inf")),
            (-float("inf"), "<=", x2, "<=", float("inf")),
            (-2, "<=", sys.u[0], "<=", 1),
            (300, "==", eco_opt.t_f(x1)),
            (0, "==", eco_opt.t_f(x2)),
        ],
    )
    return ocp, sys


def solve(control_intervals=30, polynomial_degree=5, tf_guess=30.0):
    ocp, sys = build_ocp()
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        tf_guess=tf_guess,
    )
    return sol, sys


# Analytic reference values (see module docstring derivation).
ANALYTIC_TF = 30.0
ANALYTIC_SWITCH_TIME = 20.0


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("tf:", sol.tf, "(analytic:", ANALYTIC_TF, ")")
