"""Van der Pol oscillator: CasADi / PROPT-style Lagrange regulation.

Problem (Andersson, Gillis, Diehl; CasADi direct-collocation demo)::

    minimize   int_0^T (x1^2 + x2^2 + u^2) dt
    subject to dx1/dt = (1 - x2^2) x1 - x2 + u
               dx2/dt = x1
               x1(0) = 0,  x2(0) = 1
               -1 <= u(t) <= 1
               T = 10  (fixed)

``x1`` is the oscillator "velocity" in the first-order form used by
CasADi (``ẋ1 = (1-x2²)x1 - x2 + u``, ``ẋ2 = x1``), ``x2`` the
position-like state, ``u`` an additive force. The running cost drives
the oscillation toward the origin without a terminal-state constraint.

No closed-form optimum. Pontryagin stationarity of the unconstrained
Hamiltonian gives ``u = -λ1/2``, projected onto ``[-1, 1]``. The
objective is regression-tracked against eco_opt's own converged value
at the default mesh; a coarse/fine mesh check guards the number.

See ``README.md`` for the PMP, the relation to the PROPT Mayer form
(auxiliary integral state, tf=5, tighter u-bounds), and what the tests
check.
"""

import eco_opt
from eco_opt import Namespace

T_HORIZON = 10.0
X1_0 = 0.0
X2_0 = 1.0
U_MIN, U_MAX = -1.0, 1.0


def vdp_model(t, x, u):
    x1, x2 = x[0], x[1]
    force = u[0]
    y = Namespace(x1=x1, x2=x2, u=force)
    return [(1.0 - x2 ** 2) * x1 - x2 + force, x1], y


def build_ocp(T=T_HORIZON):
    sys = eco_opt.System(states=2, controls=1, model=vdp_model, name="vdp")
    ocp = eco_opt.Ocp()
    ocp.min([sys.y.x1 ** 2 + sys.y.x2 ** 2 + sys.y.u ** 2])
    ocp.st(
        system=sys,
        constraints=[
            (U_MIN, "<=", sys.u[0], "<=", U_MAX),
            (0.0, "==", eco_opt.t_0(sys.t)),
            (T, "==", eco_opt.t_f(sys.t)),
            (X1_0, "==", eco_opt.t_0(sys.y.x1)),
            (X2_0, "==", eco_opt.t_0(sys.y.x2)),
        ],
    )
    return ocp, sys


def solve(T=T_HORIZON, control_intervals=40, polynomial_degree=4):
    ocp, sys = build_ocp(T=T)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: X1_0, sys.x[1]: X2_0, sys.u[0]: 0.0},
        ipopt_options={"max_iter": 2000, "tol": 1e-8},
        tf_guess=T,
    )
    return sol, sys


# Default-mesh Lagrange cost (control_intervals=40, degree=4, Legendre).
# Frozen so a silent transcription change fails the test.
REGRESSION_OBJECTIVE = 2.8878976747752554


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("J =", sol.objective)
    print("u min/max", float(sol.u[:, 0].min()), float(sol.u[:, 0].max()))
    print("x(0)", sol.x[0], "x(tf)", sol.x[-1])
