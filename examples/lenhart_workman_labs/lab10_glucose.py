"""Lenhart and Workman, Lab 10 (glucose–insulin).

The book’s T=20 case is outside the GTT model’s validity; default here
is the physically reasonable T=0.2 (book 16.2).
"""

import casadi as cs

import eco_opt
from eco_opt import Namespace

# Book (16.2): T = 0.2 days ≈ 4.8 h, inside GTT range.
DEFAULT = dict(a=1.0, b=1.0, c=1.0, x10=0.75, A=2.0, l=0.5, T=0.2)


def _model(a, b, c):
    def f(t, x, u):
        x1, x2 = x[0], x[1]
        return cs.vertcat(-a * x1 - b * x2, -c * x2 + u[0]), Namespace(
            x1=x1, x2=x2, u=u[0]
        )

    return f


def build_ocp(**kw):
    p = dict(DEFAULT)
    p.update({k: float(v) for k, v in kw.items() if v is not None})
    sys = eco_opt.System(
        states=2, controls=1, model=_model(p["a"], p["b"], p["c"]), name="lab10"
    )
    ocp = eco_opt.Ocp()
    ocp.min([p["A"] * (sys.y.x1 - p["l"]) ** 2 + sys.y.u**2])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (p["T"], "==", eco_opt.t_f(sys.t)),
            (p["x10"], "==", eco_opt.t_0(sys.y.x1)),
            (0.0, "==", eco_opt.t_0(sys.y.x2)),
        ],
    )
    return ocp, sys, p


def solve(control_intervals=40, polynomial_degree=3, **kw):
    ocp, sys, p = build_ocp(**kw)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: p["x10"], sys.x[1]: 0.0, sys.u[0]: 0.0},
        tf_guess=p["T"],
    )
    return sol, sys, p


if __name__ == "__main__":
    sol, _, p = solve()
    print(sol)
    print("x1[0], x1[-1], target", float(sol.x[0, 0]), float(sol.x[-1, 0]), p["l"])
