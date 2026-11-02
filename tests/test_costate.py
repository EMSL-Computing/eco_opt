"""Shadow prices (costates) and Hamiltonian recovered from the NLP multipliers.

Checked against closed-form costates (microbial growth, Manzoni et al. 2023;
one- and two-state linear-quadratic problems; a Mayer-only problem) and
against finite-difference sensitivities dJ*/dx for a two-state model, both
at t = 0 and, by re-solving the remaining horizon (Bellman's principle), at
interior times.
"""
import os
import sys as _sys

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import eco_opt
from eco_opt import Namespace

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(example, module):
    path = os.path.join(HERE, "examples", example)
    if path not in _sys.path:
        _sys.path.insert(0, path)
    return __import__(module)


@pytest.mark.parametrize("points", ["legendre", "radau"])
def test_microbial_costate_matches_pontryagin(points):
    mg = _load("microbial_growth", "microbial_growth")
    ocp, _ = mg.build_ocp()
    sol = ocp.solve(
        control_intervals=200,
        polynomial_degree=5,
        collocation_points=points,
        tf_guess=2.0,
        ipopt_options={"max_iter": 5000},
    )
    assert sol.success
    assert sol.costate.shape == (200 * 5, 1)
    assert sol.costate_time.shape == (200 * 5,)

    # lam(t) = lam_T exp(gamma (t - T)), with lam_T from H(T) = 0.
    lam_an = mg.analytical_costate(sol.costate_time, sol.tf)
    assert np.max(np.abs(sol.costate[:, 0] - lam_an)) < 1e-3 * lam_an[-1]

    # Free final time and autonomous dynamics: H = 0 along the optimum.
    assert np.max(np.abs(sol.hamiltonian)) < 1e-3 * sol.objective


def _lq_model(t, x, u):
    return u[0], Namespace(x=x[0], u=u[0])


@pytest.mark.parametrize("points", ["legendre", "radau"])
def test_lq_min_costate_and_constant_hamiltonian(points):
    s = eco_opt.System(states=1, controls=1, model=_lq_model)
    ocp = eco_opt.Ocp()
    ocp.min([0.5 * (s.y.x ** 2 + s.y.u ** 2)])
    ocp.st(
        system=s,
        constraints=[
            (0, "==", eco_opt.t_0(s.t)),
            (1, "==", eco_opt.t_f(s.t)),
            (1, "==", eco_opt.t_0(s.y.x)),
        ],
    )
    sol = ocp.solve(control_intervals=50, polynomial_degree=4, collocation_points=points)
    assert sol.success
    lam_an = np.sinh(1.0 - sol.costate_time) / np.cosh(1.0)
    assert np.max(np.abs(sol.costate[:, 0] - lam_an)) < 1e-4
    # H = (x^2 - lam^2)/2 is constant; at t = 1, lam = 0 and x = 1/cosh(1).
    assert np.max(np.abs(sol.hamiltonian - 0.5 / np.cosh(1.0) ** 2)) < 1e-4


def _double_integrator(t, x, u):
    return [x[1], u[0]], Namespace(x1=x[0], x2=x[1], u=u[0])


@pytest.mark.parametrize("points", ["legendre", "radau"])
def test_two_state_lq_costate_matches_riccati(points):
    # min int_0^1 (x1^2 + u^2)/2, x1' = x2, x2' = u, free final state:
    # lam(t) = P(t) x(t), with -P' = A'P + PA - PBB'P + Q and P(1) = 0.
    s = eco_opt.System(states=2, controls=1, model=_double_integrator)
    ocp = eco_opt.Ocp()
    ocp.min([0.5 * (s.y.x1 ** 2 + s.y.u ** 2)])
    ocp.st(
        system=s,
        constraints=[
            (0, "==", eco_opt.t_0(s.t)),
            (1, "==", eco_opt.t_f(s.t)),
            (1, "==", eco_opt.t_0(s.y.x1)),
            (0, "==", eco_opt.t_0(s.y.x2)),
        ],
    )
    sol = ocp.solve(control_intervals=50, polynomial_degree=4, collocation_points=points)
    assert sol.success

    A = np.array([[0.0, 1.0], [0.0, 0.0]])
    B = np.array([[0.0], [1.0]])
    Q = np.diag([1.0, 0.0])

    def riccati(t, p):
        P = p.reshape(2, 2)
        return -(A.T @ P + P @ A - P @ B @ B.T @ P + Q).ravel()

    def closed_loop(t, x):
        return (A - B @ B.T @ P_of(t)) @ x

    ric = solve_ivp(riccati, (1.0, 0.0), np.zeros(4), dense_output=True, rtol=1e-11, atol=1e-13)
    P_of = lambda t: ric.sol(t).reshape(2, 2)
    traj = solve_ivp(closed_loop, (0.0, 1.0), [1.0, 0.0], dense_output=True,
                     rtol=1e-11, atol=1e-13)
    lam_an = np.array([P_of(t) @ traj.sol(t) for t in sol.costate_time])
    assert np.max(np.abs(sol.costate - lam_an)) < 1e-4 * np.max(np.abs(lam_an))


def test_mayer_only_costate_is_terminal_gradient():
    s = eco_opt.System(states=1, controls=1, model=_lq_model)
    ocp = eco_opt.Ocp()
    ocp.min([eco_opt.t_f(0.5 * s.y.x ** 2), 0.5 * s.y.u ** 2])
    ocp.st(
        system=s,
        constraints=[
            (0, "==", eco_opt.t_0(s.t)),
            (1, "==", eco_opt.t_f(s.t)),
            (1, "==", eco_opt.t_0(s.y.x)),
        ],
    )
    sol = ocp.solve(control_intervals=20, polynomial_degree=3)
    # u = -lam, lam = x(1) = 1 - lam  =>  lam = 1/2 everywhere.
    assert np.allclose(sol.costate[:, 0], 0.5, atol=1e-8)


def test_two_state_costate_matches_finite_difference_sensitivity():
    ld = _load("opt_lignin", "lignin_decomposition")
    ocp, sy = ld.compile_ocp(control_intervals=60, polynomial_degree=4)
    opts = {"max_iter": 3000, "tol": 1e-12}
    base = ocp.optimize(ipopt_options=opts)
    assert base.success

    CT0 = ld.CT_0_DEFAULT
    Co0 = CT0 * ld.ARC0_LOOP

    def J_at(CT, Co):
        ocp.set_constraints(
            [(CT, "==", eco_opt.t_0(sy.y.CT)), (Co, "==", eco_opt.t_0(sy.y.Co))]
        )
        return ocp.optimize(ipopt_options=opts).objective

    eps = 1e-4 * CT0
    fd = np.array(
        [
            (J_at(CT0 + eps, Co0) - J_at(CT0 - eps, Co0)) / (2 * eps),
            (J_at(CT0, Co0 + eps) - J_at(CT0, Co0 - eps)) / (2 * eps),
        ]
    )
    # Linear extrapolation from the first two collocation points to t = 0.
    (t1, t2), (l1, l2) = base.costate_time[:2], base.costate[:2]
    lam0 = l1 - (l2 - l1) / (t2 - t1) * t1
    assert np.allclose(lam0, fd, rtol=1e-3)

    # Autonomous, fixed final time: H is constant along the optimum.
    assert np.ptp(base.hamiltonian) < 1e-6


def test_two_state_costate_matches_tail_sensitivity_at_interior_times():
    # Bellman: lam(t_k) = dV/dx at x*(t_k), where V solves the remaining
    # horizon [t_k, T] (shifted to start at 0; the model is autonomous).
    ld = _load("opt_lignin", "lignin_decomposition")
    N, d = 60, 4
    ocp, _ = ld.compile_ocp(control_intervals=N, polynomial_degree=d)
    opts = {"max_iter": 3000, "tol": 1e-12}
    base = ocp.optimize(ipopt_options=opts)
    assert base.success

    for k in (15, 45):
        tk, (CT, Co) = base.time[k], base.x[k]
        tail, ty = ld.compile_ocp(
            TT=ld.TT_DEFAULT - tk, control_intervals=N - k, polynomial_degree=d
        )

        def V_at(CT, Co):
            tail.set_constraints(
                [(CT, "==", eco_opt.t_0(ty.y.CT)), (Co, "==", eco_opt.t_0(ty.y.Co))]
            )
            out = tail.optimize(ipopt_options=opts)
            assert out.success
            return out.objective

        eps = 1e-4 * CT
        fd = np.array(
            [
                (V_at(CT + eps, Co) - V_at(CT - eps, Co)) / (2 * eps),
                (V_at(CT, Co + eps) - V_at(CT, Co - eps)) / (2 * eps),
            ]
        )
        lam = [np.interp(tk, base.costate_time, base.costate[:, j]) for j in range(2)]
        assert np.allclose(lam, fd, rtol=1e-3)


def test_multiphase_results_have_no_costate():
    def m(t, x, u):
        return u[0], Namespace(x=x[0], u=u[0])

    s1 = eco_opt.System(states=1, controls=1, model=m)
    s2 = eco_opt.System(states=1, controls=1, model=m)
    p1 = eco_opt.Phase(
        system=s1,
        constraints=[
            (0, "==", eco_opt.t_0(s1.t)),
            (0.5, "==", eco_opt.t_f(s1.t)),
            (1, "==", eco_opt.t_0(s1.y.x)),
            (-1, "<=", s1.y.u, "<=", 1),
        ],
    )
    p2 = eco_opt.Phase(
        system=s2,
        constraints=[(1, "==", eco_opt.t_f(s2.t)), (-1, "<=", s2.y.u, "<=", 1)],
    )
    ocp = eco_opt.Ocp()
    ocp.min([eco_opt.t_f(s2.y.x ** 2)])
    ocp.st(phases=[p1, p2])
    sol = ocp.solve(control_intervals=10, polynomial_degree=3)
    assert sol.costate is None and sol.hamiltonian is None
