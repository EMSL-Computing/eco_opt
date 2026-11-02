"""Microbial growth optimization for a finite food source
(Manzoni, Chakrawal & Ledder 2023, https://doi.org/10.3389/fevo.2023.1094269).

Problem (free terminal time)::

    maximize   int_0^T G(u) dt
    subject to Cdot = -u - gamma*C + mu*G(u),  C(0) = 1,  0 <= C <= 1
               G(u) = emax*beta*(u - r)/(u + beta) - k*u
               0 <= u <= 10,  T free (T >= 0)

``C`` is the substrate/resource pool, ``u`` the uptake rate (control), and
``G`` the growth rate. Growth recycles a fraction ``mu`` of biomass back into
the resource pool; resources are additionally lost via first-order decay at
rate ``gamma``.

Reference
---------
For the ``k = 0`` case used here, Manzoni, Chakrawal and Ledder (2023)
give a closed-form optimum (their ``optidec5``). See
``analytical_solution``.
"""

from __future__ import annotations

import numpy as np

import eco_opt
from eco_opt import Namespace

# Parameters (k=0 baseline case)
BETA = 0.5
GAMMA = 0.5
R = 0.1
EMAX = 0.5
MU = 0.5
K = 0.0

U_MIN, U_MAX = 0.0, 10.0
C_MIN, C_MAX = 0.0, 1.0
C0 = 1.0


def growth_rate(u, emax=EMAX, beta=BETA, r=R, k=K):
    return emax * beta * (u - r) / (u + beta) - k * u


def analytical_solution(
    c0=C0,
    beta=BETA,
    gamma=GAMMA,
    emax=EMAX,
    rho=R,
    mu=MU,
    n=401,
):
    """Closed-form optimum for ``k = 0`` (Ledder's ``optidec5``).

    Nondimensionalisation follows Manzoni, Chakrawal and Ledder (2023).
    ``rho`` is the maintenance parameter ``r`` in the growth function.
    Returns time, substrate, uptake, growth, terminal time, and
    ``J = \\int G\\,dt`` (trapezoid on the returned grid). ``eta`` is
    the efficiency diagnostic from ``optidec5``, not the OCP objective.
    """
    C0n = gamma * c0 / beta
    rho_n = rho / beta
    mu_n = emax * mu
    u_tau = rho_n + np.sqrt((1.0 + rho_n) * rho_n)
    v_tau = 1.0 + u_tau
    msq = mu_n * (1.0 + rho_n)
    m = np.sqrt(msq)
    u0 = rho_n + np.sqrt((1.0 + rho_n) * (C0n + rho_n))
    T_tau = np.log((C0n + 2.0 * u0 + 1.0 - mu_n) / (2.0 * u_tau + 1.0 - mu_n))
    T = np.linspace(0.0, T_tau, n)
    S = T_tau - T
    vsq = (v_tau**2 - msq) * np.exp(S) + msq
    v = np.sqrt(vsq)
    U = v - 1.0
    G = (U - rho_n) / v
    C = vsq / (1.0 + rho_n) - 2.0 * v + 1.0
    if mu == 0:
        eta = (emax / C0n) * (
            T_tau - 2.0 * (1.0 + rho_n) * (1.0 - np.exp(-T_tau / 2.0)) / v_tau
        )
    else:
        qd = np.sqrt(v_tau**2 - m**2 * (1.0 - np.exp(-T_tau))) + m * np.exp(
            -T_tau / 2.0
        )
        eta = (emax * (1.0 - mu) / C0n) * (
            T_tau - 2.0 * (1.0 + rho_n) * np.log((v_tau + m) / qd) / m
        )
    t = T / gamma
    g = emax * beta * G
    return {
        "t": t,
        "c": beta * C / gamma,
        "u": beta * U,
        "g": g,
        "tf": T_tau / gamma,
        "J": float(np.trapezoid(g, t)),
        "eta": eta,
    }


def microbial_growth_model(t, x, u):
    Cs = x[0]
    uptake = u[0]
    G = growth_rate(uptake)
    dx = -uptake - GAMMA * Cs + MU * G
    y = Namespace(Cs=Cs, uptake=uptake, growth=G)
    return dx, y  # no free tf upper bound (mirrors commented-out t_f(Cs)==0)


def build_ocp():
    sys = eco_opt.System(states=1, controls=1, model=microbial_growth_model, name="microbe")
    time = sys.t
    microbes = sys.y

    ocp = eco_opt.Ocp()
    ocp.max([growth_rate(microbes.uptake)])  # Lagrange: maximize integral of G
    ocp.st(
        system=sys,
        constraints=[
            (C_MIN, "<=", microbes.Cs, "<=", C_MAX),
            (U_MIN, "<=", microbes.uptake, "<=", U_MAX),
            (0, "==", eco_opt.t_0(time)),
            (C0, "==", eco_opt.t_0(microbes.Cs)),
            (0, "<=", eco_opt.t_f(time), "<=", float("inf")),
        ],
    )
    return ocp, sys


def solve(control_intervals=200, polynomial_degree=5, tf_guess=2.0):
    ocp, sys = build_ocp()
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="radau",
        ipopt_options={"max_iter": 5000},
        tf_guess=tf_guess,
    )
    return sol, sys


# Regression-tracked numerical baseline; the analytical benchmark is provided
# by analytical_solution(). Computed at control_intervals=200, degree=5,
# radau collocation, default (1.0) initial guess, tf_guess=2.0.
REGRESSION_OBJECTIVE = 0.15184969799373915
REGRESSION_TF = 1.5181904140431468


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective (integral of G):", sol.objective)
    print("final time (resource exhaustion):", sol.tf)
    print("C(tf):", sol.x[-1, 0])
