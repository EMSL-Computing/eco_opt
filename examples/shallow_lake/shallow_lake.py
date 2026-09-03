"""Shallow-lake phosphorus: ecological tipping point and two management locals.

Dimensionless lake (Mäler, Xepapadeas & de Zeeuw 2003; the q=2, m=1, r=1
Hill recycling of Carpenter-type models)::

    dP/dt = L - b P + P^2 / (1 + P^2)

P is water-column phosphorus, L is external loading (the control).
Internal recycling P^2/(1+P^2) is weak when P is small and saturates when
P is large, which can give two stable equilibria (clear vs turbid) and
hysteresis in L.

Two OCPs share that f:

1. Economic management (classic): trade agricultural benefit ln L against
   damage c P^2, discounted. For some (b, c, P(0)) two locals exist —
   manage toward a cleaner lake or accept a dirtier one.
2. Restoration: start turbid, force P(T) <= P_clean, stay close to a
   desired long-run loading. Optimal L is typically cut hard, then eased
   up after the lake has crossed the threshold.
"""

import numpy as np
import casadi as cs

import eco_opt
from eco_opt import Namespace

# Physical lake (Scheffer / Carpenter dimensionless recycling).
B = 0.52

# Economic management.
RHO = 0.03
C_DAMAGE = 1.0
T_MANAGE = 80.0
P0_MANAGE = 1.0
L_MIN, L_MAX = 1e-3, 1.5

# Restoration.
T_RESTORE = 40.0
P0_RESTORE = 1.8
P_CLEAN = 0.40
L_DESIRED = 0.05
L_RESTORE_MAX = 0.40

N_MANAGE = 50
N_RESTORE = 40
DEGREE = 3


def recycling(P):
    return P**2 / (1.0 + P**2)


def loading_at_equilibrium(P, b=B):
    """L such that dP/dt = 0 at this P (constant loading)."""
    P = np.asarray(P, dtype=float)
    return b * P - recycling(P)


def equilibria_at_loading(L, b=B):
    """Real nonnegative roots of dP/dt = 0 at constant L.

    Cubic: b P^3 - (L+1) P^2 + b P - L = 0.
    """
    roots = np.roots([b, -(L + 1.0), b, -L])
    real = np.real(roots[np.isclose(np.imag(roots), 0.0)])
    return np.sort(real[real >= -1e-12])


def lake_model(t, x, u):
    P, L = x[0], u[0]
    dP = L - B * P + recycling(P)
    y = Namespace(P=P, L=L, recycling=recycling(P))
    return dP, y


def build_management_ocp(P0=P0_MANAGE, T=T_MANAGE, c=C_DAMAGE, rho=RHO):
    sys = eco_opt.System(states=1, controls=1, model=lake_model, name="lake")
    ocp = eco_opt.Ocp()
    ocp.max([cs.exp(-rho * sys.t) * (cs.log(sys.y.L) - c * sys.y.P**2)])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (T, "==", eco_opt.t_f(sys.t)),
            (P0, "==", eco_opt.t_0(sys.y.P)),
            (0.0, "<=", sys.y.P),
            (L_MIN, "<=", sys.y.L, "<=", L_MAX),
        ],
    )
    return ocp, sys


def solve_management(P0=P0_MANAGE, L_guess=0.05, T=T_MANAGE, c=C_DAMAGE,
                     control_intervals=N_MANAGE, polynomial_degree=DEGREE):
    ocp, sys = build_management_ocp(P0=P0, T=T, c=c)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: P0, sys.u[0]: L_guess},
        tf_guess=T,
        ipopt_options={"max_iter": 3000, "tol": 1e-7},
    )
    return sol, sys


def build_restoration_ocp(
    P0=P0_RESTORE, T=T_RESTORE, P_clean=P_CLEAN, L_desired=L_DESIRED
):
    sys = eco_opt.System(states=1, controls=1, model=lake_model, name="restore")
    ocp = eco_opt.Ocp()
    ocp.min([(sys.y.L - L_desired) ** 2])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (T, "==", eco_opt.t_f(sys.t)),
            (P0, "==", eco_opt.t_0(sys.y.P)),
            (0.0, "<=", sys.y.P),
            (0.0, "<=", sys.y.L, "<=", L_RESTORE_MAX),
            (P_clean, ">=", eco_opt.t_f(sys.y.P)),
        ],
    )
    return ocp, sys


def solve_restoration(
    P0=P0_RESTORE,
    T=T_RESTORE,
    P_clean=P_CLEAN,
    L_desired=L_DESIRED,
    L_guess=0.0,
    control_intervals=N_RESTORE,
    polynomial_degree=DEGREE,
):
    ocp, sys = build_restoration_ocp(
        P0=P0, T=T, P_clean=P_clean, L_desired=L_desired
    )
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: P0, sys.u[0]: L_guess},
        tf_guess=T,
        ipopt_options={"max_iter": 3000, "tol": 1e-7},
    )
    return sol, sys


# Frozen at N=50, degree 3, T=80, P0=1, c=1, b=0.52, rho=0.03.
REGRESSION_J_LOW_GUESS = -135.223
REGRESSION_J_HIGH_GUESS = -132.850


if __name__ == "__main__":
    print("equilibria at L=0.05:", equilibria_at_loading(0.05))
    sol_lo, _ = solve_management(L_guess=0.05)
    sol_hi, _ = solve_management(L_guess=0.50)
    print("management low-L guess ", sol_lo)
    print("  P(T)=", sol_lo.x[-1, 0], "J=", sol_lo.objective)
    print("management high-L guess", sol_hi)
    print("  P(T)=", sol_hi.x[-1, 0], "J=", sol_hi.objective)
    sol_r, _ = solve_restoration()
    print("restoration", sol_r)
    print("  P(T)=", sol_r.x[-1, 0], "L min/mean",
          float(sol_r.u[:, 0].min()), float(sol_r.u[:, 0].mean()))
