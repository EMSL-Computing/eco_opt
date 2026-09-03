"""Cohen (1971) terminal-yield allocation: bang-bang vegetative → reproductive.

Problem (fixed season, linear vegetative production)::

    maximize   x2(T)
    subject to dx1/dt = u * x1
               dx2/dt = (1 - u) * x1
               x1(0) = X10 > 0,  x2(0) = 0
               0 <= u(t) <= 1,   t in [0, T]

``x1`` is vegetative biomass, ``x2`` reproductive biomass, ``u`` the
fraction of photosynthate kept in the vegetative body. Instantaneous
production equals vegetative mass (relative growth rate 1 in these units).

Analytic solution
-----------------
The Hamiltonian is affine in ``u``, so the optimum is bang-bang. With
relative growth rate 1 the switching function changes sign when the
remaining season equals one unit of time:

    u*(t) = 1  for t in [0, t*)     (pure vegetative)
    u*(t) = 0  for t in [t*, T]     (pure reproductive)
    t*    = T - 1

States::

    x1(t) = X10 * exp(t)                 0 <= t <= t*
    x1(t) = X10 * exp(t*)                t* <= t <= T
    x2(t) = 0                            0 <= t <= t*
    x2(t) = X10 * exp(t*) * (t - t*)     t* <= t <= T
    x2(T) = X10 * exp(T - 1)

A useful corollary: ``x1(T) = x2(T)`` when vegetative and reproductive
conversion efficiencies are equal (both 1 here).

Default parameters: T=5, X10=1  ⇒  t*=4,  x2(T)=exp(4).

See ``README.md`` in this directory for the Pontryagin derivation, the
link to Lindh et al. (2016) exponential growth, and the King–Roughgarden
graded-allocation companion.
"""

import math

import casadi as cs
import numpy as np

import eco_opt
from eco_opt import Namespace

T_SEASON = 5.0
X10 = 1.0
ANALYTIC_SWITCH = T_SEASON - 1.0
ANALYTIC_YIELD = X10 * math.exp(T_SEASON - 1.0)


def plant_model(t, x, u):
    veg, repro = x[0], x[1]
    alloc = u[0]
    dx = cs.vertcat(alloc * veg, (1.0 - alloc) * veg)
    y = Namespace(veg=veg, repro=repro, alloc=alloc)
    return dx, y


def build_ocp(T=T_SEASON, x10=X10):
    sys = eco_opt.System(states=2, controls=1, model=plant_model, name="cohen")
    ocp = eco_opt.Ocp()
    ocp.max(eco_opt.t_f(sys.y.repro))
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "<=", sys.u[0], "<=", 1.0),
            (0.0, "<=", sys.y.veg),
            (0.0, "<=", sys.y.repro),
            (0.0, "==", eco_opt.t_0(sys.t)),
            (T, "==", eco_opt.t_f(sys.t)),
            (x10, "==", eco_opt.t_0(sys.y.veg)),
            (0.0, "==", eco_opt.t_0(sys.y.repro)),
        ],
    )
    return ocp, sys


def solve(T=T_SEASON, x10=X10, control_intervals=80, polynomial_degree=3):
    ocp, sys = build_ocp(T=T, x10=x10)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: x10, sys.x[1]: 0.0, sys.u[0]: 0.5},
        ipopt_options={"max_iter": 2000, "tol": 1e-8},
        tf_guess=T,
    )
    return sol, sys


def analytic_switch(T=T_SEASON):
    return T - 1.0


def analytic_yield(T=T_SEASON, x10=X10):
    return x10 * math.exp(T - 1.0)


def analytic_trajectories(t, T=T_SEASON, x10=X10):
    """Closed-form states and control on an arbitrary time grid."""
    t = np.asarray(t, dtype=float)
    t_star = analytic_switch(T)
    veg_at_switch = x10 * np.exp(t_star)
    veg = np.where(t <= t_star, x10 * np.exp(t), veg_at_switch)
    repro = np.where(t <= t_star, 0.0, veg_at_switch * (t - t_star))
    alloc = np.where(t < t_star, 1.0, 0.0)
    return veg, repro, alloc


def flowering_time_from_ocp(sol, thresh=0.5):
    """Last interval-start at which the control is still vegetative (u > thresh)."""
    u = np.asarray(sol.u[:, 0])
    t = np.asarray(sol.time_intervals[:-1])
    high = u > thresh
    if not np.any(high):
        return float(t[0])
    return float(t[int(np.nonzero(high)[0][-1])])


if __name__ == "__main__":
    sol, sys = solve()
    tF = flowering_time_from_ocp(sol)
    print(sol)
    print("analytic t* =", ANALYTIC_SWITCH, "OCP tF =", tF)
    print("analytic x2(T) =", ANALYTIC_YIELD, "OCP =", sol.objective)
    print("x1(T) =", sol.x[-1, 0], "x2(T) =", sol.x[-1, 1])
