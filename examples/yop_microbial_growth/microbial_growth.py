"""Microbial growth optimization for a finite food source (Manzoni et al. 2022).

Direct Python/eco_opt translation of the MATLAB Yop model in
``examples/OptLignin/YOP benchmark/Manzoni et al. 2022/{GrowthModel.m,
yop_opt.m, YOP_microbial_growth.m}`` (Chakrawal/Manzoni,
https://doi.org/10.3389/fevo.2023.1094269).

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
No closed-form objective is given in the paper; the MATLAB script instead
solves the *same* problem three ways (forward-backward sweep, ``bvp4c``
indirect shooting via Pontryagin's maximum principle, and the Yop/CasADi
direct-collocation solve reproduced here, ``indirect_bvp4c.m`` /
``indirect_forward_backward_sweep.m`` / ``YOP_microbial_growth.m``) and
cross-checks them against one another. The indirect (``bvp4c``) costate
equations derive a stationarity condition ``dH/du = 0`` on the Hamiltonian
``H = G(u) + p*(-u - gamma*C + mu*G(u))`` that is costate- (hence
time-) dependent, so the optimal control here is a smooth trajectory over
the horizon, not a constant singular arc.

We validate physically instead of against a closed-form number:
* the free terminal time must be driven by resource exhaustion, i.e.
  ``C(tf) -> 0`` (finite substrate optimal-uptake problem);
* the control must respect its box bounds ``0 <= u <= 10`` throughout;
* the objective is regression-tracked (``REGRESSION_OBJECTIVE`` below) to
  catch silent future regressions in the transcription.
"""

import casadi as cs

import eco_opt
from eco_opt.system import Namespace

# Parameters (see YOP_microbial_growth.m, k=0 baseline case)
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


# Regression-tracked value: no external closed-form objective available (see
# module docstring); computed with eco_opt at control_intervals=200, degree=5,
# radau collocation, default (1.0) initial guess, tf_guess=2.0.
REGRESSION_OBJECTIVE = 0.15184969799373915
REGRESSION_TF = 1.5181904140431468


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective (integral of G):", sol.objective)
    print("final time (resource exhaustion):", sol.tf)
    print("C(tf):", sol.x[-1, 0])
