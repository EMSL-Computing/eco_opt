"""Mixed path constraints must hold at interval-start nodes.

Radau collocation does not include tau=0. Piecewise-constant u is defined
on the left of each interval, so g(x,u) at those knots is not implied by
collocation-point inequalities.
"""

import numpy as np

import eco_opt
from eco_opt import Namespace


def _solve_u_le_x(collocation_points="radau"):
    import casadi as cs

    def model(t, x, u):
        dx = cs.vertcat(1.0)
        y = Namespace(x=x[0], u=u[0], slack=x[0] - u[0])
        return dx, y

    sys = eco_opt.System(states=1, controls=1, model=model, name="u_le_x")
    y = sys.y
    ocp = eco_opt.Ocp()
    ocp.max([y.u])
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(sys.t)),
            (1, "==", eco_opt.t_f(sys.t)),
            (0, "==", eco_opt.t_0(y.x)),
            (0, "<=", y.slack),
            (-2, "<=", y.u, "<=", 2),
            (0, "<=", y.x),
        ],
    )
    sol = ocp.solve(
        control_intervals=8,
        polynomial_degree=3,
        collocation_points=collocation_points,
        ipopt_options={"print_level": 0, "sb": "yes"},
        tf_guess=1.0,
    )
    return sol


def test_mixed_path_constraint_holds_at_t0_radau():
    sol = _solve_u_le_x("radau")
    assert sol.success, sol.return_status
    # x(0)=0 ⇒ u(0)<=0. Without interval-start path rows, U_0 can track
    # x at the first collocation point (>0).
    assert sol.u[0, 0] <= sol.x[0, 0] + 1e-6


def test_mixed_path_constraint_holds_at_every_interval_start():
    sol = _solve_u_le_x("radau")
    assert sol.success, sol.return_status
    n = sol.u.shape[0]
    slack = sol.x[:n, 0] - sol.u[:, 0]
    assert np.all(slack >= -1e-6), slack.min()


def test_mixed_path_constraint_holds_at_t0_legendre():
    sol = _solve_u_le_x("legendre")
    assert sol.success, sol.return_status
    assert sol.u[0, 0] <= sol.x[0, 0] + 1e-6
