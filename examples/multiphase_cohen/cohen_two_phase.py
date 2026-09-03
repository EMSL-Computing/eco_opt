"""Cohen (1971) as a two-phase OCP: vegetative then reproductive.

Single-phase Cohen (``examples/cohen_king_roughgarden``) lets ``u(t)``
jump inside one collocation grid. Here the bangs are *baked into the
dynamics of two phases* and the optimizer chooses only the switch time:

    Phase 1 (vegetative):   dx1/dt = x1,  dx2/dt = 0
    Phase 2 (reproductive): dx1/dt = 0,   dx2/dt = x1
    link: t_f^(1) = t_0^(2),  x^(1)(t_s) = x^(2)(t_s)
    maximize x2(T),  T fixed, t_s free

Analytic: t_s* = T - 1,  x2(T) = X10 * exp(T - 1).
Defaults T=5, X10=1  ⇒  t_s*=4,  yield=exp(4).
"""

import math

import casadi as cs
import numpy as np

import eco_opt
from eco_opt import Namespace, Phase

T_SEASON = 5.0
X10 = 1.0
ANALYTIC_SWITCH = T_SEASON - 1.0
ANALYTIC_YIELD = X10 * math.exp(T_SEASON - 1.0)


def vegetative_model(t, x, u):
    veg, repro = x[0], x[1]
    dx = cs.vertcat(veg, 0.0 * repro)
    y = Namespace(veg=veg, repro=repro)
    return dx, y


def reproductive_model(t, x, u):
    veg, repro = x[0], x[1]
    dx = cs.vertcat(0.0 * veg, veg)
    y = Namespace(veg=veg, repro=repro)
    return dx, y


def build_ocp(T=T_SEASON, x10=X10):
    sys_v = eco_opt.System(
        states=2, controls=0, model=vegetative_model, name="cohen_veg"
    )
    sys_r = eco_opt.System(
        states=2, controls=0, model=reproductive_model, name="cohen_repro"
    )
    p1 = Phase(
        sys_v,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys_v.t)),
            (x10, "==", eco_opt.t_0(sys_v.y.veg)),
            (0.0, "==", eco_opt.t_0(sys_v.y.repro)),
            (0.0, "<=", sys_v.y.veg),
            (0.0, "<=", sys_v.y.repro),
        ],
        control_intervals=20,
        name="vegetative",
    )
    p2 = Phase(
        sys_r,
        constraints=[
            (T, "==", eco_opt.t_f(sys_r.t)),
            (0.0, "<=", sys_r.y.veg),
            (0.0, "<=", sys_r.y.repro),
        ],
        control_intervals=20,
        name="reproductive",
    )
    ocp = eco_opt.Ocp()
    ocp.max(eco_opt.t_f(sys_r.y.repro))
    ocp.st(phases=[p1, p2], continuity=True)
    return ocp, (sys_v, sys_r)


def _solve_options(systems, T, x10, tf_guess):
    return dict(
        control_intervals=20,
        polynomial_degree=3,
        collocation_points="legendre",
        initial_guess={
            systems[0].x[0]: x10,
            systems[0].x[1]: 0.0,
            systems[1].x[0]: x10 * math.exp(T - 1.0),
            systems[1].x[1]: 0.0,
        },
        ipopt_options={"max_iter": 2000, "tol": 1e-8},
        tf_guess=tf_guess,
    )


def solve(T=T_SEASON, x10=X10, tf_guess=None):
    if tf_guess is None:
        tf_guess = T  # equal knots; the switch time is a decision variable
    ocp, systems = build_ocp(T=T, x10=x10)
    sol = ocp.solve(**_solve_options(systems, T, x10, tf_guess))
    return sol, systems


def optimize(T=T_SEASON, x10=X10, tf_guess=None):
    """Same problem as ``solve``, using ``build`` then ``optimize``."""
    if tf_guess is None:
        tf_guess = T
    ocp, systems = build_ocp(T=T, x10=x10)
    opts = _solve_options(systems, T, x10, tf_guess)
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


if __name__ == "__main__":
    sol, systems = optimize()
    print(sol)
    print("n_phases", len(sol.phases))
    print("switch t_s =", sol.switch_times, "analytic", ANALYTIC_SWITCH)
    print("yield", sol.objective, "analytic", ANALYTIC_YIELD)
    print("phase 1 [t0, tf]", sol.phases[0].t0, sol.phases[0].tf)
    print("phase 2 [t0, tf]", sol.phases[1].t0, sol.phases[1].tf)
    print("x2(T)", sol.x[-1, 1], "x1(T)", sol.x[-1, 0])
