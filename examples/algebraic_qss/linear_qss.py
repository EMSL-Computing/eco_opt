"""Linear QSS microbes: the same OCP as a reduced ODE and as a DAE.

    max  ∫ M dt
    Cdot = -u
    0 = (emax - r u) u C/(K+C) - d M     (growth = death)
    C(0)=1, T=10, 0≤u≤umax, C≥0, M≥0

The algebraic residual is linear in M, so M = formula(C,u) is unique.
``solve_ode`` substitutes that formula into J (what eco_opt could always
do). ``solve_dae`` keeps M as an algebraic unknown with 0=g.
"""

import casadi as cs

import eco_opt
from eco_opt import Namespace

T_HORIZON = 10.0
C0 = 1.0
EMAX = 0.5
R = 0.4
K = 0.2
D = 0.1
UMAX = 0.5


def biomass(C, u):
    growth = (EMAX - R * u) * u * C / (K + C)
    return growth / D


def ode_model(t, x, u):
    C = x[0]
    uptake = u[0]
    M = biomass(C, uptake)
    y = Namespace(C=C, M=M, u=uptake)
    return -uptake, y


def dae_model(t, x, u, z):
    C = x[0]
    uptake = u[0]
    M = z[0]
    growth = (EMAX - R * uptake) * uptake * C / (K + C)
    g = growth - D * M
    y = Namespace(C=C, M=M, u=uptake)
    return -uptake, g, y


def _constraints(sys, C, M, uptake):
    return [
        (0.0, "==", eco_opt.t_0(sys.t)),
        (T_HORIZON, "==", eco_opt.t_f(sys.t)),
        (C0, "==", eco_opt.t_0(C)),
        (0.0, "<=", C),
        (0.0, "<=", M),
        (0.0, "<=", uptake, "<=", UMAX),
    ]


def build_ode():
    sys = eco_opt.System(states=1, controls=1, model=ode_model, name="qss_ode")
    ocp = eco_opt.Ocp()
    ocp.max([sys.y.M])
    ocp.st(system=sys, constraints=_constraints(sys, sys.y.C, sys.y.M, sys.y.u))
    return ocp, sys


def build_dae():
    sys = eco_opt.System(
        states=1, controls=1, algebraics=1, model=dae_model, name="qss_dae"
    )
    ocp = eco_opt.Ocp()
    ocp.max([sys.y.M])
    ocp.st(system=sys, constraints=_constraints(sys, sys.y.C, sys.y.M, sys.y.u))
    return ocp, sys


def _solve(ocp, sys, is_dae=False):
    guess = {sys.x[0]: C0, sys.u[0]: 0.25}
    if is_dae:
        guess[sys.z[0]] = float(biomass(C0, 0.25))
    return ocp.solve(
        control_intervals=24,
        polynomial_degree=3,
        collocation_points="radau",
        initial_guess=guess,
        tf_guess=T_HORIZON,
        ipopt_options={"max_iter": 2000, "tol": 1e-8},
    )


def solve_ode():
    ocp, sys = build_ode()
    return _solve(ocp, sys, is_dae=False), sys


def solve_dae():
    ocp, sys = build_dae()
    return _solve(ocp, sys, is_dae=True), sys


if __name__ == "__main__":
    sol_ode, _ = solve_ode()
    sol_dae, _ = solve_dae()
    print("ODE", sol_ode)
    print("DAE", sol_dae)
    print("J ode/dae", sol_ode.objective, sol_dae.objective)
