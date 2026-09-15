"""Lenhart and Workman, Lab 8 (HIV reverse-transcriptase therapy), (14.1).

u=0 is maximum therapy, u=1 is no therapy (infectivity scales with u).

    maximize   int_0^T [ A T - (1-u)^2 ] dt
    T'  = s/(1+V) - m1 T + r T (1 - (T+Ti)/Tmax) - u k V T
    Ti' = u k V T - m2 Ti
    V'  = N m2 Ti - m3 V
    0 <= u <= 1
"""

import casadi as cs

import eco_opt
from eco_opt import Namespace

DEFAULT = dict(
    s=10.0,
    m1=0.02,
    m2=0.5,
    m3=4.4,
    r=0.03,
    Tmax=1500.0,
    k=0.000024,
    N=300.0,
    T0=800.0,
    Ti0=0.04,
    V0=1.5,
    A=0.05,
    Tfinal=20.0,
)


def _model(s, m1, m2, m3, r, Tmax, k, Npart):
    def f(t, x, u):
        T, Ti, V = x[0], x[1], x[2]
        ctrl = u[0]
        dT = (
            s / (1.0 + V)
            - m1 * T
            + r * T * (1.0 - (T + Ti) / Tmax)
            - ctrl * k * V * T
        )
        dTi = ctrl * k * V * T - m2 * Ti
        dV = Npart * m2 * Ti - m3 * V
        return [dT, dTi, dV], Namespace(T=T, Ti=Ti, V=V, u=ctrl)

    return f


def build_ocp(**kw):
    p = dict(DEFAULT)
    p.update({k: float(v) for k, v in kw.items() if v is not None})
    sys = eco_opt.System(
        states=3,
        controls=1,
        model=_model(p["s"], p["m1"], p["m2"], p["m3"], p["r"], p["Tmax"], p["k"], p["N"]),
        name="lab8",
    )
    ocp = eco_opt.Ocp()
    ocp.max([p["A"] * sys.y.T - (1.0 - sys.y.u) ** 2])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (p["Tfinal"], "==", eco_opt.t_f(sys.t)),
            (p["T0"], "==", eco_opt.t_0(sys.y.T)),
            (p["Ti0"], "==", eco_opt.t_0(sys.y.Ti)),
            (p["V0"], "==", eco_opt.t_0(sys.y.V)),
            (0.0, "<=", sys.y.T, "<=", cs.inf),
            (0.0, "<=", sys.y.Ti, "<=", cs.inf),
            (0.0, "<=", sys.y.V, "<=", cs.inf),
            (0.0, "<=", sys.y.u, "<=", 1.0),
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
            sys.x[0]: p["T0"],
            sys.x[1]: p["Ti0"],
            sys.x[2]: p["V0"],
            sys.u[0]: 0.5,
        },
        tf_guess=p["Tfinal"],
    )
    return sol, sys, p


if __name__ == "__main__":
    sol, _, _ = solve()
    print(sol)
    print("u[0], u[-1]", float(sol.u[0, 0]), float(sol.u[-1, 0]))
