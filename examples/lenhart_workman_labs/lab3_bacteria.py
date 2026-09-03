"""Lenhart and Workman, Lab 3 (bacteria and nutrient).

Problem::

    maximize   C x(1) - int_0^1 u^2 dt
    subject to xdot = r x + A u x - B u^2 exp(-x),  x(0) = x0 > 0
               A, B, C >= 0

Default book values (7.1): r=1, A=1, B=12, C=1, x0=1.

Analytic facts from the book:
- A = 0 implies u* ≡ 0 (nutrient has no benefit).
- Terminal payoff: the adjoint satisfies λ(1) = C, not λ(1) = 0.
No closed-form trajectory.
"""

import casadi as cs

import eco_opt
from eco_opt import Namespace

T_HORIZON = 1.0
DEFAULT = dict(r=1.0, A=1.0, B=12.0, C=1.0, x0=1.0)


def _model(r, A, B):
    def bacteria_model(t, x, u):
        state, ctrl = x[0], u[0]
        dx = r * state + A * ctrl * state - B * ctrl**2 * cs.exp(-state)
        y = Namespace(x=state, u=ctrl)
        return dx, y

    return bacteria_model


def build_ocp(r=None, A=None, B=None, C=None, x0=None):
    p = dict(DEFAULT)
    if r is not None:
        p["r"] = float(r)
    if A is not None:
        p["A"] = float(A)
    if B is not None:
        p["B"] = float(B)
    if C is not None:
        p["C"] = float(C)
    if x0 is not None:
        p["x0"] = float(x0)

    sys = eco_opt.System(
        states=1,
        controls=1,
        model=_model(p["r"], p["A"], p["B"]),
        name="lab3",
    )
    ocp = eco_opt.Ocp()
    ocp.max([eco_opt.t_f(p["C"] * sys.y.x), -sys.y.u**2])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (T_HORIZON, "==", eco_opt.t_f(sys.t)),
            (p["x0"], "==", eco_opt.t_0(sys.y.x)),
            (0.0, "<=", sys.y.x, "<=", cs.inf),
            (-cs.inf, "<=", sys.y.u, "<=", cs.inf),
        ],
    )
    return ocp, sys, p


def solve(
    r=None,
    A=None,
    B=None,
    C=None,
    x0=None,
    control_intervals=40,
    polynomial_degree=3,
):
    ocp, sys, p = build_ocp(r=r, A=A, B=B, C=C, x0=x0)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: p["x0"], sys.u[0]: 0.0},
        tf_guess=T_HORIZON,
    )
    return sol, sys, p


if __name__ == "__main__":
    sol, _sys, p = solve()
    print(sol)
    print("params", p)
    print("u[0], u[-1]:", float(sol.u[0, 0]), float(sol.u[-1, 0]))
    print("x[0], x[-1]:", float(sol.x[0, 0]), float(sol.x[-1, 0]))
