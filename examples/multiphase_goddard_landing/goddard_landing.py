"""Goddard rocket landing: two-phase max-ascent then soft landing.

Yop example ``goddardLanding`` (yoptimization.com). Same dimensional
rocket as ``examples/goddard_benchmark/``, transcribed as two
sequential phases:

    Phase 1 (ascent):  maximize h(t_s),  t_s free
    Phase 2 (landing): h(T)=0, v(T)=0,   T free
    link: time and (v, h, m) continuous at t_s
    0 <= F <= Fm on both phases;  mf <= m <= m0;  h >= 0

The Yop listing maximises height at the end of phase 1 and uses a dummy
``min t_f(0)`` on phase 2 (no second objective). Soft landing is
enforced by the terminal constraints only. The prose problem statement
also writes ``m(T)=mf``; the listed Yop code does not pin terminal mass,
only the path bound ``mf <= m <= m0``. This file follows the code.

No closed-form altitude. Peak height is regression-tracked after a
verified solve.
"""

import casadi as cs
import numpy as np

import eco_opt
from eco_opt import Namespace, Phase

M0 = 214.839
MF = 67.9833
FM = 9.525515
D0 = 0.01227
BETA = 0.145e-3
C_EXH = 2060.0
G0 = 9.81
R0 = 6.371e6


def goddard_model(t, x, u):
    v, h, m = x[0], x[1], x[2]
    F = u[0]
    D = D0 * cs.exp(-BETA * h)
    F_D = cs.sign(v) * D * v ** 2
    g = G0 * (R0 / (R0 + h)) ** 2
    dv = (F * C_EXH - F_D) / m - g
    dh = v
    dm = -F
    dx = [dv, dh, dm]
    y = Namespace()
    y.rocket = Namespace(velocity=v, height=h, mass=m, fuelMassFlow=F)
    y.drag = Namespace(coefficient=D, force=F_D)
    y.gravity = g
    return dx, y


def build_ocp():
    sys1 = eco_opt.System(
        states=3, controls=1, model=goddard_model, name="goddard_ascent"
    )
    sys2 = eco_opt.System(
        states=3, controls=1, model=goddard_model, name="goddard_landing"
    )
    r1, r2 = sys1.y.rocket, sys2.y.rocket

    p1 = Phase(
        sys1,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys1.t)),
            (0.0, "==", eco_opt.t_0(r1.velocity)),
            (0.0, "==", eco_opt.t_0(r1.height)),
            (M0, "==", eco_opt.t_0(r1.mass)),
            (0.0, "<=", r1.height),
            (MF, "<=", r1.mass, "<=", M0),
            (0.0, "<=", r1.fuelMassFlow, "<=", FM),
        ],
        name="ascent",
    )
    p2 = Phase(
        sys2,
        constraints=[
            (0.0, "<=", eco_opt.t_0(sys2.t)),
            (0.0, "<=", r2.height),
            (MF, "<=", r2.mass, "<=", M0),
            (0.0, "<=", r2.fuelMassFlow, "<=", FM),
            (0.0, "==", eco_opt.t_f(r2.velocity)),
            (0.0, "==", eco_opt.t_f(r2.height)),
        ],
        name="landing",
    )
    ocp = eco_opt.Ocp()
    # Height at the end of phase 1 (Yop: ocp(1).max(t_f(height))).
    ocp.max(eco_opt.t_f(r1.height))
    ocp.st(phases=[p1, p2], continuity=True)
    return ocp, (sys1, sys2)


def _solve_options(systems, control_intervals, polynomial_degree, tf_guess, ipopt_options):
    sys1, sys2 = systems
    opts = {"max_iter": 4000, "tol": 1e-8, "nlp_scaling_method": "gradient-based"}
    if ipopt_options:
        opts.update(ipopt_options)
    return dict(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="radau",
        initial_guess={
            sys1.x[0]: 200.0,
            sys1.x[1]: 5000.0,
            sys1.x[2]: 0.5 * (M0 + MF),
            sys1.u[0]: 0.5 * FM,
            sys2.x[0]: -100.0,
            sys2.x[1]: 5000.0,
            sys2.x[2]: 0.5 * (M0 + MF),
            sys2.u[0]: 0.5 * FM,
        },
        ipopt_options=opts,
        tf_guess=tf_guess,
    )


def solve(
    control_intervals=50,
    polynomial_degree=3,
    tf_guess=None,
    ipopt_options=None,
):
    if tf_guess is None:
        # Ascent then descent; overall horizon free. Equal knots from 0.
        tf_guess = [0.0, 40.0, 80.0]
    ocp, systems = build_ocp()
    sol = ocp.solve(
        **_solve_options(
            systems, control_intervals, polynomial_degree, tf_guess, ipopt_options
        )
    )
    return sol, systems


def optimize(
    control_intervals=50,
    polynomial_degree=3,
    tf_guess=None,
    ipopt_options=None,
):
    """Same problem as ``solve``, using ``build`` then ``optimize``."""
    if tf_guess is None:
        tf_guess = [0.0, 40.0, 80.0]
    ocp, systems = build_ocp()
    opts = _solve_options(
        systems, control_intervals, polynomial_degree, tf_guess, ipopt_options
    )
    ocp.build(
        control_intervals=opts["control_intervals"],
        polynomial_degree=opts["polynomial_degree"],
        collocation_points=opts["collocation_points"],
        initial_guess=opts["initial_guess"],
        tf_guess=opts["tf_guess"],
    )
    sol = ocp.optimize(
        initial_guess=opts["initial_guess"],
        ipopt_options=opts["ipopt_options"],
        tf_guess=opts["tf_guess"],
        warm_start=False,
    )
    return sol, systems


# Peak altitude at end of phase 1, and overall tf, at the default mesh
# (control_intervals=50, polynomial_degree=3, radau).
REGRESSION_PEAK_HEIGHT = 82083.12418041468
REGRESSION_TF = 289.25163273395356


if __name__ == "__main__":
    sol, systems = optimize()
    print(sol)
    print("n_phases", len(sol.phases))
    print("switch t_s", sol.switch_times, "tf", sol.tf)
    print("peak h (phase-1 terminal / objective)", sol.objective)
    print("h(T), v(T), m(T)", sol.x[-1])
    print("phase-1 terminal v,h,m", sol.phases[0].x[-1])
    print("min h", float(np.min(sol.x[:, 1])), "max h", float(np.max(sol.x[:, 1])))
    print("F min/max", float(sol.u[:, 0].min()), float(sol.u[:, 0].max()))
