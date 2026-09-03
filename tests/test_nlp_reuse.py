"""Build once, retune bounds/parameters, re-optimize (no new transcription)."""
import casadi as cs
import numpy as np
import pytest

import eco_opt
from eco_opt.system import Namespace


def _energy_ocp(x_final=1.0):
    def model(t, x, u):
        p, v = x[0], x[1]
        a = u[0]
        dx = cs.vertcat(v, a)
        y = Namespace(position=p, velocity=v, accel=a)
        return dx, y

    sys = eco_opt.System(states=2, controls=1, model=model, name="di")
    y = sys.y
    ocp = eco_opt.Ocp()
    ocp.min([y.accel ** 2])
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(y.position)),
            (0, "==", eco_opt.t_0(y.velocity)),
            (x_final, "==", eco_opt.t_f(y.position)),
            (0, "==", eco_opt.t_f(y.velocity)),
            (0.0, "==", eco_opt.t_0(sys.t)),
            (1.0, "==", eco_opt.t_f(sys.t)),
            (-10, "<=", y.accel, "<=", 10),
        ],
    )
    return ocp, sys


def test_build_optimize_matches_solve():
    ocp, sys = _energy_ocp(1.0)
    sol_one = ocp.solve(control_intervals=16, polynomial_degree=3, tf_guess=1.0)

    ocp2, sys2 = _energy_ocp(1.0)
    ocp2.build(control_intervals=16, polynomial_degree=3, tf_guess=1.0)
    sol_two = ocp2.optimize(warm_start=False, tf_guess=1.0)

    assert sol_one.success and sol_two.success
    assert abs(sol_one.objective - sol_two.objective) < 1e-8
    assert abs(sol_two.x[-1, 0] - 1.0) < 1e-5


def test_set_constraints_retunes_terminal_without_rebuild():
    ocp, sys = _energy_ocp(1.0)
    ocp.build(control_intervals=16, polynomial_degree=3, tf_guess=1.0)
    sol1 = ocp.optimize(warm_start=False, tf_guess=1.0)
    solver_id = id(ocp._compiled.solver)

    ocp.set_constraints([(2.0, "==", eco_opt.t_f(sys.y.position))])
    sol2 = ocp.optimize()

    assert sol1.success and sol2.success
    assert id(ocp._compiled.solver) == solver_id
    assert abs(sol1.x[-1, 0] - 1.0) < 1e-4
    assert abs(sol2.x[-1, 0] - 2.0) < 1e-4
    assert sol2.objective > sol1.objective


def test_path_box_retune_bryson_style():
    def model(t, x, u):
        dx = cs.vertcat(x[1], u[0])
        y = Namespace(position=x[0], speed=x[1], acceleration=u[0])
        return dx, y

    sys = eco_opt.System(states=2, controls=1, model=model, name="trolley")
    y = sys.y
    ocp = eco_opt.Ocp()
    ocp.min([0.5 * y.acceleration ** 2])
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(y.position)),
            (1, "==", eco_opt.t_0(y.speed)),
            (1, "==", eco_opt.t_f(sys.t)),
            (0, "==", eco_opt.t_f(y.position)),
            (-1, "==", eco_opt.t_f(y.speed)),
            (1.0 / 9.0, ">=", y.position),
        ],
    )
    ocp.build(control_intervals=20, polynomial_degree=5, tf_guess=1.0)
    sol_loose = ocp.optimize(warm_start=False, tf_guess=1.0)
    assert sol_loose.success
    assert abs(sol_loose.objective - 4.0) / 4.0 < 5e-3

    ocp.set_constraints([(1.0 / 20.0, ">=", y.position)])
    sol_tight = ocp.optimize()
    assert sol_tight.success
    assert sol_tight.x[:, 0].max() <= 1.0 / 20.0 + 2e-3
    assert sol_tight.objective > sol_loose.objective


def test_parameter_pin_sweep_reuses_nlp():
    def model(t, x, u, p):
        dx = -p[0] * x[0] + u[0]
        y = Namespace(x=x[0], u=u[0], p=p[0])
        return dx, y

    sys = eco_opt.System(states=1, controls=1, parameters=1, model=model, name="decay")
    ocp = eco_opt.Ocp()
    ocp.min([sys.y.u ** 2])
    ocp.st(
        system=sys,
        constraints=[
            (1.0, "==", eco_opt.t_0(sys.y.x)),
            (0.2, "==", eco_opt.t_f(sys.y.x)),
            (0.0, "==", eco_opt.t_0(sys.t)),
            (1.0, "==", eco_opt.t_f(sys.t)),
            (-5, "<=", sys.y.u, "<=", 5),
            (1.0, "==", sys.p[0]),
        ],
    )
    ocp.build(control_intervals=12, polynomial_degree=3, tf_guess=1.0)
    sol1 = ocp.optimize(warm_start=False, tf_guess=1.0)
    solver_id = id(ocp._compiled.solver)

    ocp.set_constraints([(3.0, "==", sys.p[0])])
    sol2 = ocp.optimize()

    assert sol1.success and sol2.success
    assert id(ocp._compiled.solver) == solver_id
    assert abs(sol1.p[0] - 1.0) < 1e-8
    assert abs(sol2.p[0] - 3.0) < 1e-8
    assert abs(sol1.objective - sol2.objective) > 1e-4


def test_new_path_expression_requires_rebuild():
    ocp, sys = _energy_ocp(1.0)
    ocp.build(control_intervals=8, polynomial_degree=3, tf_guess=1.0)
    with pytest.raises(ValueError, match="not in the compiled NLP"):
        ocp.set_constraints([(0.0, "<=", sys.y.position ** 2)])


def test_multiphase_set_constraints_requires_rebuild():
    def model(t, x, u):
        return cs.vertcat(x[0])

    sys_a = eco_opt.System(states=1, controls=0, model=model, name="a")
    sys_b = eco_opt.System(states=1, controls=0, model=model, name="b")
    p1 = eco_opt.Phase(
        sys_a,
        constraints=[
            (0, "==", eco_opt.t_0(sys_a.t)),
            (1, "==", eco_opt.t_0(sys_a.x[0])),
        ],
        control_intervals=4,
    )
    p2 = eco_opt.Phase(
        sys_b,
        constraints=[(1, "==", eco_opt.t_f(sys_b.t))],
        control_intervals=4,
    )
    ocp = eco_opt.Ocp()
    ocp.min(eco_opt.t_f(sys_b.x[0]))
    ocp.st(phases=[p1, p2])
    ocp.build(control_intervals=4, polynomial_degree=3, tf_guess=1.0)
    with pytest.raises(NotImplementedError, match="single-phase"):
        ocp.set_constraints([(0.5, "==", eco_opt.t_0(sys_a.x[0]))])
