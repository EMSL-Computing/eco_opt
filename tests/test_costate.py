"""Shadow prices (costates) and Hamiltonian recovered from the NLP multipliers.

Checked against closed-form costates (microbial growth, Manzoni et al. 2023;
a linear-quadratic problem; a Mayer-only problem) and against a
finite-difference sensitivity dJ*/dx0 for a two-state model.
"""
import os
import sys as _sys

import numpy as np
import pytest

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
    a = mg.EMAX * mg.BETA
    u_T = mg.R + np.sqrt(mg.R ** 2 + mg.R * mg.BETA)
    g_u = a * (mg.R + mg.BETA) / (u_T + mg.BETA) ** 2
    lam_T = g_u / (1.0 - mg.MU * g_u)
    lam_an = lam_T * np.exp(mg.GAMMA * (sol.costate_time - sol.tf))
    assert np.max(np.abs(sol.costate[:, 0] - lam_an)) < 1e-3 * lam_T

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
