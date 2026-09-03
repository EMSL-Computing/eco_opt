"""Lenhart and Workman, Lab 13 (predator–prey pesticide, isoperimetric).

    minimize   N1(T) + (A/2) int_0^T u^2 dt
    subject to Lotka–Volterra with pesticide, 0 <= u <= M,
               int_0^T u dt = B  (extra state z, z(T)=B)
"""

import casadi as cs

import eco_opt
from eco_opt import Namespace

DEFAULT = dict(
    d1=0.1, d2=0.01, N10=5.0, N20=2.0, M=1.0, A=1.0, B=1.0, T=5.0
)


def _model(d1, d2):
    def f(t, x, u):
        N1, N2, z = x[0], x[1], x[2]
        ctrl = u[0]
        dN1 = (1.0 - N2) * N1 - d1 * N1 * ctrl
        dN2 = (N1 - 1.0) * N2 - d2 * N2 * ctrl
        dz = ctrl
        return cs.vertcat(dN1, dN2, dz), Namespace(N1=N1, N2=N2, z=z, u=ctrl)

    return f


def build_ocp(**kw):
    p = dict(DEFAULT)
    p.update({k: float(v) for k, v in kw.items() if v is not None})
    if p["B"] > p["M"] * p["T"] + 1e-12:
        raise ValueError("need B <= M T")
    sys = eco_opt.System(
        states=3, controls=1, model=_model(p["d1"], p["d2"]), name="lab13"
    )
    ocp = eco_opt.Ocp()
    ocp.min([eco_opt.t_f(sys.y.N1), 0.5 * p["A"] * sys.y.u**2])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (p["T"], "==", eco_opt.t_f(sys.t)),
            (p["N10"], "==", eco_opt.t_0(sys.y.N1)),
            (p["N20"], "==", eco_opt.t_0(sys.y.N2)),
            (0.0, "==", eco_opt.t_0(sys.y.z)),
            (p["B"], "==", eco_opt.t_f(sys.y.z)),
            (0.0, "<=", sys.y.N1, "<=", cs.inf),
            (0.0, "<=", sys.y.N2, "<=", cs.inf),
            (0.0, "<=", sys.y.u, "<=", p["M"]),
        ],
    )
    return ocp, sys, p


def solve(control_intervals=50, polynomial_degree=3, **kw):
    ocp, sys, p = build_ocp(**kw)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={
            sys.x[0]: p["N10"],
            sys.x[1]: p["N20"],
            sys.x[2]: 0.0,
            sys.u[0]: p["B"] / p["T"],
        },
        tf_guess=p["T"],
    )
    return sol, sys, p


if __name__ == "__main__":
    sol, _, p = solve()
    print(sol)
    print("z(T) vs B", float(sol.x[-1, 2]), p["B"])
