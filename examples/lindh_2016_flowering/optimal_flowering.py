"""Lindh et al. (2016) "Constrained growth flips the direction of optimal
phenological responses among annual plants" (New Phytologist 209:1591-1599,
doi:10.1111/nph.13706).

Reproduces the paper's central theoretical result (Theorem 1, and the
key qualitative Fig. 2 pattern): the optimal control for allocating
biomass production between vegetative growth V(t) and reproductive
growth R(t),

    dV/dt = u(t) F(V),      V(t0) = V0
    dR/dt = (1 - u(t)) F(V), R(t0) = 0
    maximize R(E)   (E = season end),  u(t) in [0, 1]

is **bang-bang**: u=1 (pure vegetative growth) until a single switching
("flowering") time tF, then u=0 (pure reproduction) until season end.
This holds for all three growth models considered (F(V) = P * U(V)):

    exponential (unconstrained): U(V) = V
    logistic    (constrained):   U(V) = V (1 - V/Vmax)
    WBE         (constrained):   U(V) = V^(3/4) (1 - (V/Vmax)^(1/4))

and, crucially, the *direction* of d(tF*)/dP flips: monotonically
increasing for exponential growth, but hump-shaped (interior maximum)
for logistic/WBE growth -- i.e. constrained-growth plants flower
*earlier* as productivity keeps rising past some threshold, unlike
unconstrained-growth plants which always flower later.

We validate this two independent ways:
  1. Solve the *continuous-relaxation* OCP (u(t) in [0,1], no bang-bang
     assumption) with eco_opt/CasADi direct collocation, and confirm
     the optimizer's u(t) numerically saturates at the box bounds with
     a single switch (Theorem 1, verified rather than assumed).
  2. Cross-check tF*(P) against a fast closed-form/1-D approach: since
     bang-bang implies W(tF) = F(V(tF)) (E - tF) with V(t) obeying the
     *un-controlled* ODE dV/dt = F(V) up to tF (analytically
     integrable for all three models), tF* = argmax_tF W(tF) can be
     found via 1-D scalar optimization -- far cheaper than solving a
     full OCP at every P, letting us sweep P densely.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

import eco_opt
from eco_opt.system import Namespace

# ---------------------------------------------------------------------------
# Growth models: F(V) = P * U(V)
# ---------------------------------------------------------------------------
VMAX = 1.0
V0 = 0.02
SEASON_START = 0.0
SEASON_END = 1.0  # T = 1


def U_exponential(V):
    return V


def U_logistic(V, Vmax=VMAX):
    return V * (1 - V / Vmax)


def U_wbe(V, Vmax=VMAX):
    return V**0.75 * (1 - (V / Vmax) ** 0.25)


GROWTH_MODELS = {
    "exponential": U_exponential,
    "logistic": U_logistic,
    "wbe": U_wbe,
}


# ---------------------------------------------------------------------------
# 1. Full OCP (bang-bang verification): relaxed control u(t) in [0,1]
# ---------------------------------------------------------------------------
def build_ocp(growth, P, t0=SEASON_START, tf=SEASON_END, v0=V0, vmax=VMAX):
    """2-state (V, R) OCP: maximize R(tf), dV=u*F(V), dR=(1-u)*F(V)."""
    U_fn = GROWTH_MODELS[growth]

    def model(t, x, u):
        V, R = x[0], x[1]
        import casadi as cs
        Vc = cs.fmax(V, 1e-9)  # guard against V=0 in fractional powers (WBE)
        F = P * U_fn(Vc, vmax) if growth != "exponential" else P * U_fn(Vc)
        dV = u[0] * F
        dR = (1 - u[0]) * F
        y = Namespace(V=V, R=R, F=F, growth_rate=F / Vc)
        return cs.vertcat(dV, dR), y

    sys = eco_opt.System(states=2, controls=1, model=model, name=f"lindh_{growth}")
    ocp = eco_opt.Ocp()
    ocp.max([eco_opt.t_f(sys.y.R)])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "<=", sys.u[0], "<=", 1.0),
            (0.0, "==", eco_opt.t_0(sys.t)),
            (tf, "==", eco_opt.t_f(sys.t)),
            (v0, "==", eco_opt.t_0(sys.y.V)),
            (0.0, "==", eco_opt.t_0(sys.y.R)),
            (0.0, "<=", sys.y.V),
        ],
    )
    return ocp, sys


def solve_ocp(growth, P, control_intervals=80, polynomial_degree=3, tf_guess=SEASON_END):
    ocp, sys = build_ocp(growth, P)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: V0, sys.x[1]: 0.0, sys.u[0]: 0.5},
        ipopt_options={"max_iter": 1000, "tol": 1e-8},
        tf_guess=tf_guess,
    )
    return sol, sys


def flowering_time_from_ocp(sol, tol=0.1):
    """Estimate the bang-bang switch time from the solved control trajectory:
    the last time index at which u is still close to 1 (before dropping to 0)."""
    u = sol.u[:, 0]
    t = sol.time_intervals[:-1]
    high = u > (1 - tol)
    if not np.any(high):
        return t[0]
    idx = np.where(high)[0][-1]
    return t[idx]


def is_bang_bang(sol, tol=0.05):
    """Check that u(t) stays within `tol` of {0,1} for (almost) the whole
    horizon, i.e. is bang-bang rather than a smoothly varying singular arc."""
    u = sol.u[:, 0]
    near_bounds = np.minimum(np.abs(u - 0.0), np.abs(u - 1.0)) < tol
    return np.mean(near_bounds) > 0.9  # allow one transition interval


# ---------------------------------------------------------------------------
# 2. Fast closed-form/1-D cross-check for tF*(P): given bang-bang control,
#    V(t) obeys the uncontrolled ODE dV/dt = F(V) for t in [0, tF], and
#    W(tF) = F(V(tF)) * (E - tF).
# ---------------------------------------------------------------------------
def integrate_V(growth, P, tF, v0=V0, vmax=VMAX):
    U_fn = GROWTH_MODELS[growth]

    def rhs(t, V):
        Vc = max(V[0], 1e-9)
        return [P * U_fn(Vc, vmax) if growth != "exponential" else P * U_fn(Vc)]

    sol = solve_ivp(rhs, [0, tF], [v0], method="RK45", rtol=1e-9, atol=1e-12, dense_output=False)
    return float(sol.y[0, -1])


def W_of_tF(growth, P, tF, E=SEASON_END, v0=V0, vmax=VMAX):
    if tF <= 0:
        VtF = v0
    else:
        VtF = integrate_V(growth, P, tF, v0, vmax)
    U_fn = GROWTH_MODELS[growth]
    Vc = max(VtF, 1e-9)
    F = P * U_fn(Vc, vmax) if growth != "exponential" else P * U_fn(Vc)
    return F * (E - tF)


def optimal_flowering_time_1d(growth, P, E=SEASON_END, v0=V0, vmax=VMAX):
    """argmax_{tF in [0,E]} W(tF) via bounded scalar optimization."""
    res = minimize_scalar(
        lambda tF: -W_of_tF(growth, P, tF, E, v0, vmax),
        bounds=(1e-6, E - 1e-6),
        method="bounded",
        options={"xatol": 1e-6},
    )
    return res.x


if __name__ == "__main__":
    for growth in GROWTH_MODELS:
        P = 2.0
        sol, sys = solve_ocp(growth, P)
        tF_ocp = flowering_time_from_ocp(sol)
        tF_1d = optimal_flowering_time_1d(growth, P)
        print(f"{growth:12s} P={P}: OCP success={sol.success} "
              f"bang-bang={is_bang_bang(sol)} tF_ocp={tF_ocp:.4f} tF_1d={tF_1d:.4f} "
              f"objective(OCP)={sol.objective:.4f}")
