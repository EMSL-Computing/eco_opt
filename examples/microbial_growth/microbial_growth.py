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

An independent Pontryagin trajectory (same necessary conditions as
Manzoni, Chakrawal and Ledder 2023) is ``analytical_solution``. The
collocation ``J`` and ``T`` are also regression-tracked.
"""

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


def analytical_solution(ds=1e-4):
    """Independent PMP trajectory for the baseline ``k=0`` parameters.

    Integrates the costate and state backward from resource exhaustion.
    Used to compare collocation against the paper's necessary conditions;
    not part of the first-solve paste.
    """
    a = EMAX * BETA
    ut = R + np.sqrt(R ** 2 + R * BETA)
    g_u = a * (R + BETA) / (ut + BETA) ** 2
    lam_t = g_u / (1.0 - MU * g_u)

    def growth(u):
        return a * (u - R) / (u + BETA) - K * u

    def u_from_lam(lam):
        return -BETA + np.sqrt(a * (R + BETA) * (MU + 1.0 / lam))

    s = 0.0
    C = 0.0
    rows = []
    for _ in range(200000):
        lam = lam_t * np.exp(-GAMMA * s)
        u = float(u_from_lam(lam))
        g = float(growth(u))
        C = C + ds * (u + GAMMA * C - MU * g)
        s += ds
        rows.append((s, C, u, g))
        if C >= C0:
            break
    arr = np.asarray(rows)
    T = float(arr[-1, 0])
    t = T - arr[::-1, 0]
    C = arr[::-1, 1]
    u = arr[::-1, 2]
    g = arr[::-1, 3]
    trapz = getattr(np, "trapezoid", None)
    if trapz is None:
        trapz = np.trapz
    J = float(trapz(g, t))
    return t, C, u, g, T, J


# Collocation at control_intervals=200, degree=5, radau, tf_guess=2.0.
REGRESSION_OBJECTIVE = 0.15184969799373915
REGRESSION_TF = 1.5181904140431468


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective (integral of G):", sol.objective)
    print("final time (resource exhaustion):", sol.tf)
    print("C(tf):", sol.x[-1, 0])
    from plot_microbial_growth import write_vs_analytical_figure

    write_vs_analytical_figure(sol)
