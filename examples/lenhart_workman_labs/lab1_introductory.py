"""Lenhart and Workman, Lab 1 (introductory scalar example).

Problem::

    maximize   int_0^1 (A x - B u^2) dt
    subject to xdot = -0.5 x^2 + C u,  x(0) = x0 > -2
               A >= 0, B > 0

Default book values: A=1, B=1, C=4, x0=1.

Analytic facts from the book (no closed-form trajectory):
- A = 0 implies u* ≡ 0 (the state is ignored).
- Trajectories depend on A/B, not A and B separately; (A,B)=(1,2)
  and (2,4) give the same (x, u).
- B = 0 and A > 0 has no optimal control (Exercise 5.1); not solved here.
"""

import casadi as cs

import eco_opt
from eco_opt import Namespace

T_HORIZON = 1.0
DEFAULT = dict(A=1.0, B=1.0, C=4.0, x0=1.0)


def _model(C):
    def lab1_model(t, x, u):
        state, ctrl = x[0], u[0]
        dx = -0.5 * state**2 + C * ctrl
        y = Namespace(x=state, u=ctrl)
        return dx, y

    return lab1_model


def build_ocp(A=None, B=None, C=None, x0=None):
    A = DEFAULT["A"] if A is None else float(A)
    B = DEFAULT["B"] if B is None else float(B)
    C = DEFAULT["C"] if C is None else float(C)
    x0 = DEFAULT["x0"] if x0 is None else float(x0)
    if B <= 0:
        raise ValueError("B must be positive (no optimum exists at B=0)")
    if A < 0:
        raise ValueError("A must be non-negative")
    if x0 <= -2:
        raise ValueError("x0 must be greater than -2")

    sys = eco_opt.System(states=1, controls=1, model=_model(C), name="lab1")
    ocp = eco_opt.Ocp()
    ocp.max([A * sys.y.x - B * sys.y.u**2])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (T_HORIZON, "==", eco_opt.t_f(sys.t)),
            (x0, "==", eco_opt.t_0(sys.y.x)),
            (-cs.inf, "<=", sys.y.x, "<=", cs.inf),
            (-cs.inf, "<=", sys.y.u, "<=", cs.inf),
        ],
    )
    return ocp, sys


def solve(
    A=None,
    B=None,
    C=None,
    x0=None,
    control_intervals=40,
    polynomial_degree=3,
):
    ocp, sys = build_ocp(A=A, B=B, C=C, x0=x0)
    x0_val = DEFAULT["x0"] if x0 is None else float(x0)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: x0_val, sys.u[0]: 0.0},
        tf_guess=T_HORIZON,
    )
    return sol, sys


if __name__ == "__main__":
    sol, _sys = solve()
    print(sol)
    print("u[0], u[-1]:", float(sol.u[0, 0]), float(sol.u[-1, 0]))
    print("x[0], x[-1]:", float(sol.x[0, 0]), float(sol.x[-1, 0]))
