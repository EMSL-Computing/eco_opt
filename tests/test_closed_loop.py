"""Open-loop replay and shrinking-horizon closed loop."""
import numpy as np
import pytest

import eco_opt
from eco_opt import Namespace


def _lin_model(t, x, u):
    dx = -x[0] + u[0]
    return dx, Namespace(x=x[0], u=u[0])


def _vdp_ocp(T=5.0):
    def model(t, x, u):
        x1, x2 = x[0], x[1]
        force = u[0]
        return [(1.0 - x2**2) * x1 - x2 + force, x1], Namespace(
            x1=x1, x2=x2, u=force
        )

    sys = eco_opt.System(states=2, controls=1, model=model, name="vdp")
    ocp = eco_opt.Ocp()
    ocp.min([sys.y.x1**2 + sys.y.x2**2 + sys.y.u**2])
    ocp.st(
        system=sys,
        constraints=[
            (-1.0, "<=", sys.u[0], "<=", 1.0),
            (0.0, "==", eco_opt.t_0(sys.t)),
            (T, "==", eco_opt.t_f(sys.t)),
            (0.0, "==", eco_opt.t_0(sys.y.x1)),
            (1.0, "==", eco_opt.t_0(sys.y.x2)),
        ],
    )
    return ocp, sys


def test_simulate_linear_ode():
    sys = eco_opt.System(states=1, controls=1, model=_lin_model, name="lin")
    t = np.linspace(0.0, 2.0, 21)
    sim = eco_opt.simulate(sys, x0=[0.0], t_grid=t, u=1.0)
    exact = 1.0 - np.exp(-t)
    assert np.max(np.abs(sim.x[:, 0] - exact)) < 1e-4


def test_replay_matches_collocation_on_vdp():
    T = 5.0
    ocp, sys = _vdp_ocp(T)
    sol = ocp.solve(
        control_intervals=20,
        polynomial_degree=3,
        collocation_points="radau",
        initial_guess={sys.x[0]: 0.0, sys.x[1]: 1.0, sys.u[0]: 0.0},
        tf_guess=T,
        ipopt_options={"print_level": 0, "max_iter": 2000},
    )
    assert sol.success
    sim = ocp.simulate(sol)
    assert sim.x.shape[0] == sol.x.shape[0]
    # Collocation is not CVODES; they should still track.
    assert np.max(np.abs(sim.x - sol.x)) < 0.05


def test_shrinking_closed_loop_matches_open_loop_vdp():
    T = 5.0
    N = 10
    dt = T / N
    ocp, sys = _vdp_ocp(T)
    sol = ocp.solve(
        control_intervals=N,
        polynomial_degree=3,
        collocation_points="radau",
        initial_guess={sys.x[0]: 0.0, sys.x[1]: 1.0, sys.u[0]: 0.0},
        tf_guess=T,
        ipopt_options={"print_level": 0, "max_iter": 2000},
    )
    assert sol.success

    ocp2, sys2 = _vdp_ocp(T)
    cl = ocp2.closed_loop(
        x0=[0.0, 1.0],
        dt=dt,
        t_final=T,
        horizon="shrinking",
        polynomial_degree=3,
        collocation_points="radau",
        ipopt_options={"print_level": 0, "max_iter": 2000},
        initial_guess={sys2.u[0]: 0.0},
    )
    assert cl.success
    assert cl.horizon == "shrinking"
    assert abs(cl.tf - T) < 1e-9
    assert cl.u.shape[0] == N
    # Same plant as the OCP: prefixes reprint the open-loop tape.
    assert np.max(np.abs(cl.u[:, 0] - sol.u[:, 0])) < 0.08
    assert np.max(np.abs(cl.x - sol.x)) < 0.08


def test_closed_loop_differs_when_plant_is_not_the_model():
    T = 5.0
    N = 10
    dt = T / N
    ocp, sys = _vdp_ocp(T)
    sol = ocp.solve(
        control_intervals=N,
        polynomial_degree=3,
        collocation_points="radau",
        initial_guess={sys.x[0]: 0.0, sys.x[1]: 1.0, sys.u[0]: 0.0},
        tf_guess=T,
        ipopt_options={"print_level": 0, "max_iter": 2000},
    )
    assert sol.success

    def kicked(t, x, u):
        x1, x2 = x[0], x[1]
        force = u[0]
        return [(1.0 - x2**2) * x1 - x2 + force + 0.35, x1], Namespace(
            x1=x1, x2=x2, u=force
        )

    plant = eco_opt.System(states=2, controls=1, model=kicked, name="vdp_kick")
    ocp2, sys2 = _vdp_ocp(T)
    cl = ocp2.closed_loop(
        x0=[0.0, 1.0],
        dt=dt,
        t_final=T,
        horizon="shrinking",
        plant=plant,
        polynomial_degree=3,
        collocation_points="radau",
        ipopt_options={"print_level": 0, "max_iter": 2000},
        initial_guess={sys2.u[0]: 0.0},
    )
    assert cl.success
    # After a few samples the applied u is not the undisturbed tape.
    assert np.max(np.abs(cl.u[:, 0] - sol.u[:, 0])) > 0.05


def test_closed_loop_rejects_multiphase():
    def model(t, x, u):
        return x[0]

    sys_a = eco_opt.System(states=1, controls=0, model=model, name="a")
    sys_b = eco_opt.System(states=1, controls=0, model=model, name="b")
    p1 = eco_opt.Phase(
        sys_a,
        constraints=[
            (0, "==", eco_opt.t_0(sys_a.t)),
            (1, "==", eco_opt.t_0(sys_a.x[0])),
        ],
    )
    p2 = eco_opt.Phase(
        sys_b,
        constraints=[(1, "==", eco_opt.t_f(sys_b.t))],
    )
    ocp = eco_opt.Ocp()
    ocp.min(eco_opt.t_f(sys_b.x[0]))
    ocp.st(phases=[p1, p2])
    with pytest.raises(NotImplementedError, match="single-phase"):
        ocp.closed_loop(x0=[1.0], dt=0.1, t_final=1.0)


def test_public_api_exports_simulation():
    assert hasattr(eco_opt, "simulate")
    assert hasattr(eco_opt, "closed_loop")
    assert hasattr(eco_opt, "SimulationResults")
    assert hasattr(eco_opt, "ClosedLoopResults")
