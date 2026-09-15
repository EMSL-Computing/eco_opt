"""Lenhart and Workman, Lab 7 (SEIR vaccination), book case (13.1).

    minimize   int_0^T (A I + u^2) dt
    S' = b N - d S - c S I - u S
    E' = c S I - (e + d) E
    I' = e E - (g + a + d) I
    N' = (b - d) N - a I
    0 <= u <= 0.9
    N(0) = S0+E0+I0+R0

R is omitted (book: other variables do not depend on R).
Incidence is mass-action cSI, not cSI/N.
"""

import casadi as cs

import eco_opt
from eco_opt import Namespace

DEFAULT = dict(
    b=0.525,
    d=0.5,
    c=0.0001,
    e=0.5,
    g=0.1,
    a=0.2,
    S0=1000.0,
    E0=100.0,
    I0=50.0,
    R0=15.0,
    A=0.1,
    T=20.0,
)


def _model(b, d, c, e, g, a):
    def f(t, x, u):
        S, E, I, N = x[0], x[1], x[2], x[3]
        v = u[0]
        dS = b * N - d * S - c * S * I - v * S
        dE = c * S * I - (e + d) * E
        dI = e * E - (g + a + d) * I
        dN = (b - d) * N - a * I
        y = Namespace(S=S, E=E, I=I, N=N, u=v)
        return [dS, dE, dI, dN], y

    return f


def build_ocp(**kw):
    p = dict(DEFAULT)
    p.update({k: float(v) for k, v in kw.items() if v is not None})
    N0 = p["S0"] + p["E0"] + p["I0"] + p["R0"]
    sys = eco_opt.System(
        states=4,
        controls=1,
        model=_model(p["b"], p["d"], p["c"], p["e"], p["g"], p["a"]),
        name="lab7",
    )
    ocp = eco_opt.Ocp()
    ocp.min([p["A"] * sys.y.I + sys.y.u**2])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (p["T"], "==", eco_opt.t_f(sys.t)),
            (p["S0"], "==", eco_opt.t_0(sys.y.S)),
            (p["E0"], "==", eco_opt.t_0(sys.y.E)),
            (p["I0"], "==", eco_opt.t_0(sys.y.I)),
            (N0, "==", eco_opt.t_0(sys.y.N)),
            (0.0, "<=", sys.y.S, "<=", cs.inf),
            (0.0, "<=", sys.y.E, "<=", cs.inf),
            (0.0, "<=", sys.y.I, "<=", cs.inf),
            (0.0, "<=", sys.y.N, "<=", cs.inf),
            (0.0, "<=", sys.y.u, "<=", 0.9),
        ],
    )
    return ocp, sys, p


def solve(control_intervals=40, polynomial_degree=3, **kw):
    ocp, sys, p = build_ocp(**kw)
    N0 = p["S0"] + p["E0"] + p["I0"] + p["R0"]
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={
            sys.x[0]: p["S0"],
            sys.x[1]: p["E0"],
            sys.x[2]: p["I0"],
            sys.x[3]: N0,
            sys.u[0]: 0.1,
        },
        tf_guess=p["T"],
    )
    return sol, sys, p


if __name__ == "__main__":
    sol, _, p = solve()
    print(sol)
    print("I[0], I[-1]", float(sol.x[0, 2]), float(sol.x[-1, 2]))
