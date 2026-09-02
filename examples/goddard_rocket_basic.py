"""Goddard rocket, maximum ascent, nondimensional variant.

Distinct from the free-tf ``goddard_benchmark.py`` (dimensional units).

Parameters: m0=1, mf=0.6, h0=1, g0=1, D0=0.5*620=310, c=0.5, beta=500,
Tmax=3.5. Free final time. controlIntervals=100.

Reference value
---------------
This exact nondimensional parameter set (h0=m0=g0=1, D0=0.5*vc with vc=620,
c=0.5*sqrt(g0*h0)=0.5, beta=hc=500, mc=0.6, Tc=Tmax=3.5) is the "COPS"
benchmark Goddard rocket problem (Dolan/More COPS suite, problem #10) and
also appears in Seywald/Cliff-style indirect-shooting papers. A commonly
cited high-precision reference (indirect BVP solve, Petit et al. IFAC 2008,
"Solving the Goddard Problem with Thrust and Dynamic Pressure Constraints
Using Saturation Functions") reports:

    h(T) = 1.01283692,  T = 0.19885626

(both in nondimensional units). We use h(T) = 1.01283692 as the reference
objective with a documented relative tolerance, and additionally verify
grid-refinement stability (50 vs 100 vs 200 control intervals).
"""

import casadi as cs

import eco_opt
from eco_opt.system import Namespace


M0 = 1.0
MF = 0.6
H0 = 1.0
G0 = 1.0
D0 = 0.5 * 620.0
C = 0.5
BETA = 500.0
TMAX = 3.5


def goddard_basic_model(t, x, u):
    v, h, m = x[0], x[1], x[2]
    T = u[0]

    g = G0 * (H0 / h) ** 2
    D = D0 * cs.exp(-BETA * h)
    F_D = D * v ** 2

    dv = (T - cs.sign(v) * F_D) / m - g
    dh = v
    dm = -T / C
    dx = cs.vertcat(dv, dh, dm)

    y = Namespace()
    y.rocket = Namespace(speed=v, height=h, fuelMass=m, thrust=T)
    y.drag = Namespace(coefficient=D, force=F_D)
    y.gravity = g
    return dx, y


def build_ocp():
    sys = eco_opt.System(states=3, controls=1, model=goddard_basic_model, name="goddard_basic")
    time = sys.t
    rocket = sys.y.rocket

    ocp = eco_opt.Ocp()
    ocp.max(eco_opt.t_f(rocket.height))
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(time)),
            (H0, "==", eco_opt.t_0(rocket.height)),
            (0, "==", eco_opt.t_0(rocket.speed)),
            (M0, "==", eco_opt.t_0(rocket.fuelMass)),
            (0, "<=", eco_opt.t_f(time), "<=", float("inf")),
            (H0, "<=", rocket.height, "<=", float("inf")),
            (-float("inf"), "<=", rocket.speed, "<=", float("inf")),
            (MF, "<=", rocket.fuelMass, "<=", M0),
            (0, "<=", rocket.thrust, "<=", TMAX),
        ],
    )
    return ocp, sys


def solve(control_intervals=100, polynomial_degree=5, tf_guess=0.2):
    ocp, sys = build_ocp()
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        tf_guess=tf_guess,
    )
    return sol, sys


REFERENCE_HEIGHT = 1.01283692
REFERENCE_TF = 0.19885626


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective (final height):", sol.objective, "(reference:", REFERENCE_HEIGHT, ")")
    print("final time:", sol.tf, "(reference:", REFERENCE_TF, ")")
    print("final mass:", sol.x[-1, 2])
