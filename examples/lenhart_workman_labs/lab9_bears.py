"""Lenhart and Workman, Lab 9 (bear metapopulation), book case (15.1).

    min int (O + cp up^2 + cf uf^2) dt
    P' = r P - (r/K) P^2 + (mf r/K)(1 - P/K) F^2 - up P
    F' = r F - (r/K) F^2 + (mp r/K)(1 - F/K) P^2 - uf F
    O' = r(1-mp) P^2/K + r(1-mf) F^2/K
         + (mf r/K^2) P F^2 + (mp r/K^2) P^2 F
    0 <= up, uf <= 1
Harvest is only on P and F, not on O.
"""

import casadi as cs

import eco_opt
from eco_opt import Namespace

DEFAULT = dict(
    r=0.1,
    K=0.75,
    mp=0.5,
    mf=0.5,
    P0=0.4,
    F0=0.2,
    O0=0.0,
    cp=10000.0,
    cf=10.0,
    T=25.0,
)


def _model(r, K, mp, mf):
    def f(t, x, u):
        P, F, O = x[0], x[1], x[2]
        up, uf = u[0], u[1]
        dP = r * P - (r / K) * P**2 + (mf * r / K) * (1.0 - P / K) * F**2 - up * P
        dF = r * F - (r / K) * F**2 + (mp * r / K) * (1.0 - F / K) * P**2 - uf * F
        dO = (
            r * (1.0 - mp) * P**2 / K
            + r * (1.0 - mf) * F**2 / K
            + (mf * r / K**2) * P * F**2
            + (mp * r / K**2) * P**2 * F
        )
        return cs.vertcat(dP, dF, dO), Namespace(P=P, F=F, O=O, up=up, uf=uf)

    return f


def build_ocp(**kw):
    p = dict(DEFAULT)
    p.update({k: float(v) for k, v in kw.items() if v is not None})
    sys = eco_opt.System(
        states=3,
        controls=2,
        model=_model(p["r"], p["K"], p["mp"], p["mf"]),
        name="lab9",
    )
    ocp = eco_opt.Ocp()
    ocp.min([sys.y.O + p["cp"] * sys.y.up**2 + p["cf"] * sys.y.uf**2])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (p["T"], "==", eco_opt.t_f(sys.t)),
            (p["P0"], "==", eco_opt.t_0(sys.y.P)),
            (p["F0"], "==", eco_opt.t_0(sys.y.F)),
            (p["O0"], "==", eco_opt.t_0(sys.y.O)),
            (0.0, "<=", sys.y.P, "<=", cs.inf),
            (0.0, "<=", sys.y.F, "<=", cs.inf),
            (0.0, "<=", sys.y.O, "<=", cs.inf),
            (0.0, "<=", sys.y.up, "<=", 1.0),
            (0.0, "<=", sys.y.uf, "<=", 1.0),
        ],
    )
    return ocp, sys, p


def solve(control_intervals=40, polynomial_degree=3, **kw):
    ocp, sys, p = build_ocp(**kw)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={
            sys.x[0]: p["P0"],
            sys.x[1]: p["F0"],
            sys.x[2]: p["O0"],
            sys.u[0]: 0.0,
            sys.u[1]: 0.1,
        },
        tf_guess=p["T"],
    )
    return sol, sys, p


if __name__ == "__main__":
    sol, _, _ = solve()
    print(sol)
    print("up max, uf max", float(sol.u[:, 0].max()), float(sol.u[:, 1].max()))
