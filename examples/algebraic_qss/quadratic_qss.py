"""QSS microbes with density-dependent death (algebraics matter).

    max  ∫ M dt
    Cdot = -u
    0 = (emax - r u) u C/(K+C) - d M - c M^2
    C(0)=1, T=10, 0≤u≤umax, C≥0, M≥0

g=0 is quadratic in M. The physical root is the positive branch of

    M = [-d + sqrt(d^2 + 4 c growth)] / (2 c)

Putting that sqrt in f is easy to get wrong (sign of the root, growth<0).
The DAE keeps 0=g and lets the NLP pick M≥0.
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
C_QUAD = 0.5
UMAX = 0.5


def growth_rate(C, u):
    return (EMAX - R * u) * u * C / (K + C)


def biomass_positive_root(C, u):
    """Closed-form QSS biomass (the root you must not mix up)."""
    g = growth_rate(C, u)
    disc = D * D + 4.0 * C_QUAD * g
    return (-D + cs.sqrt(disc)) / (2.0 * C_QUAD)


def ode_model(t, x, u):
    C = x[0]
    uptake = u[0]
    M = biomass_positive_root(C, uptake)
    y = Namespace(C=C, M=M, u=uptake)
    return -uptake, y


def dae_model(t, x, u, z):
    C = x[0]
    uptake = u[0]
    M = z[0]
    g = growth_rate(C, uptake) - D * M - C_QUAD * M ** 2
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
    sys = eco_opt.System(states=1, controls=1, model=ode_model, name="qssq_ode")
    ocp = eco_opt.Ocp()
    ocp.max([sys.y.M])
    ocp.st(system=sys, constraints=_constraints(sys, sys.y.C, sys.y.M, sys.y.u))
    return ocp, sys


def build_dae():
    sys = eco_opt.System(
        states=1, controls=1, algebraics=1, model=dae_model, name="qssq_dae"
    )
    ocp = eco_opt.Ocp()
    ocp.max([sys.y.M])
    ocp.st(system=sys, constraints=_constraints(sys, sys.y.C, sys.y.M, sys.y.u))
    return ocp, sys


def _solve(ocp, sys, is_dae=False):
    M0 = float(biomass_positive_root(C0, 0.25))
    guess = {sys.x[0]: C0, sys.u[0]: 0.25}
    if is_dae:
        guess[sys.z[0]] = M0
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
    sol_dae, sys = solve_dae()
    print("ODE", sol_ode)
    print("DAE", sol_dae)
    print("J ode/dae", sol_ode.objective, sol_dae.objective)
    print("M(0) dae", sol_dae.z[0, 0], "formula", float(biomass_positive_root(C0, sol_dae.u[0, 0])))
