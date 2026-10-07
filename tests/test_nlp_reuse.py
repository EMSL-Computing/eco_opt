"""Build once, retune bounds/parameters, re-optimize (no new transcription)."""
import numpy as np
import pytest

import eco_opt
from eco_opt.constraints import classify_constraint
from eco_opt import Namespace


def _energy_ocp(x_final=1.0):
    def model(t, x, u):
        p, v = x[0], x[1]
        a = u[0]
        y = Namespace(position=p, velocity=v, accel=a)
        return [v, a], y

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
        y = Namespace(position=x[0], speed=x[1], acceleration=u[0])
        return [x[1], u[0]], y

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
    sol_tight = ocp.optimize(warm_start=False, tf_guess=1.0)
    assert sol_tight.success
    assert sol_tight.x[:, 0].max() <= 1.0 / 20.0 + 2e-3
    assert sol_tight.objective > sol_loose.objective

    ocp.set_constraints([(1.0 / 5.0, ">=", y.position)])
    sol_loose_again = ocp.optimize(warm_start=False, tf_guess=1.0)
    assert sol_loose_again.success
    assert sol_loose_again.x[:, 0].max() <= 1.0 / 5.0 + 2e-3
    assert sol_loose_again.objective < sol_tight.objective


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


def test_set_constraints_updates_ocp_constraints_and_survives_rebuild():
    ocp, sys = _energy_ocp(1.0)
    ocp.build(control_intervals=12, polynomial_degree=3, tf_guess=1.0)
    ocp.optimize(warm_start=False, tf_guess=1.0)

    ocp.set_constraints([(2.0, "==", eco_opt.t_f(sys.y.position))])

    finals = [
        classify_constraint(sys, c)
        for c in ocp.constraints
        if classify_constraint(sys, c).kind == "box_final"
        and classify_constraint(sys, c).slot == ("x", 0)
    ]
    assert len(finals) == 1
    assert finals[0].lb == 2.0 and finals[0].ub == 2.0

    accel = [
        classify_constraint(sys, c)
        for c in ocp.constraints
        if classify_constraint(sys, c).kind == "box"
        and classify_constraint(sys, c).slot == ("u", 0)
    ]
    assert len(accel) == 1
    assert accel[0].lb == -10.0 and accel[0].ub == 10.0

    ocp.build(control_intervals=20, polynomial_degree=3, tf_guess=1.0)
    sol = ocp.optimize(warm_start=False, tf_guess=1.0)
    assert sol.success
    assert abs(sol.x[-1, 0] - 2.0) < 1e-4
    assert abs(sol.u).max() <= 10.0 + 1e-6


def test_path_bound_survives_rebuild():
    def model(t, x, u):
        y = Namespace(position=x[0], speed=x[1], acceleration=u[0])
        return [x[1], u[0]], y

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
    ocp.build(control_intervals=16, polynomial_degree=3, tf_guess=1.0)
    ocp.optimize(warm_start=False, tf_guess=1.0)
    ocp.set_constraints([(1.0 / 20.0, ">=", y.position)])
    ocp.build(control_intervals=24, polynomial_degree=3, tf_guess=1.0)
    sol = ocp.optimize(warm_start=False, tf_guess=1.0)
    assert sol.success
    assert sol.x[:, 0].max() <= 1.0 / 20.0 + 3e-3


def test_new_path_expression_requires_rebuild():
    ocp, sys = _energy_ocp(1.0)
    ocp.build(control_intervals=8, polynomial_degree=3, tf_guess=1.0)
    with pytest.raises(ValueError, match="not in the compiled NLP"):
        ocp.set_constraints([(0.0, "<=", sys.y.position ** 2)])


def test_multiphase_set_constraints_requires_rebuild():
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


# -- nlpsol_options pass-through ------------------------------------------
def _solve_di(**kw):
    ocp, _ = _energy_ocp(1.0)
    ocp.build(control_intervals=16, polynomial_degree=3, tf_guess=1.0)
    return ocp, ocp.optimize(warm_start=False, **kw)


def test_nlpsol_options_default_unchanged():
    print("eco_opt imported from", eco_opt.__file__)
    ocp, sol_a = _solve_di()
    _, sol_b = _solve_di(nlpsol_options=None)
    assert float(sol_a.objective) == float(sol_b.objective)
    assert np.array_equal(np.asarray(sol_a.x), np.asarray(sol_b.x))
    from eco_opt.ocp import _solver_opts
    assert _solver_opts(None) == {"ipopt": {"print_level": 0}, "print_time": 0}


def test_nlpsol_options_expand_matches_default():
    _, sol_a = _solve_di()
    ocp, sol_b = _solve_di(nlpsol_options={"expand": True})
    assert ocp._compiled.solver.stats()["success"]
    assert abs(float(sol_a.objective) - float(sol_b.objective)) < 1e-8
    assert np.max(np.abs(np.asarray(sol_a.x) - np.asarray(sol_b.x))) < 1e-8


def test_nlpsol_options_ipopt_merge_and_rebuild():
    from eco_opt.ocp import _solver_opts
    o = _solver_opts({"max_iter": 5, "tol": 1e-6},
                     {"expand": True, "ipopt": {"max_iter": 99, "sb": "yes"}})
    assert o["expand"] is True and o["print_time"] == 0
    assert o["ipopt"] == {"max_iter": 5, "tol": 1e-6, "sb": "yes", "print_level": 0}
    # nlpsol print_time / print_level overrides are honoured
    o = _solver_opts(None, {"print_time": 1, "ipopt": {"print_level": 3}})
    assert o["print_time"] == 1 and o["ipopt"]["print_level"] == 3
    # inputs not mutated
    nl = {"ipopt": {"tol": 1e-7}}
    _solver_opts({"max_iter": 3}, nl)
    assert nl == {"ipopt": {"tol": 1e-7}}

    # both routes reach IPOPT: max_iter=1 via either route stops early
    for kw in ({"ipopt_options": {"max_iter": 1}},
               {"nlpsol_options": {"ipopt": {"max_iter": 1}}}):
        ocp, _ = _solve_di(**kw)
        assert ocp._compiled.solver.stats()["iter_count"] <= 1
    # ipopt_options wins over nlpsol_options["ipopt"]
    ocp, _ = _solve_di(ipopt_options={"max_iter": 1},
                       nlpsol_options={"ipopt": {"max_iter": 1000}})
    assert ocp._compiled.solver.stats()["iter_count"] <= 1

    # solver is rebuilt when nlpsol_options change between calls
    ocp, _ = _solve_di()
    s1 = ocp._compiled.solver
    ocp.optimize(nlpsol_options={"expand": True})
    assert ocp._compiled.solver is not s1
    s2 = ocp._compiled.solver
    ocp.optimize(nlpsol_options={"expand": True})
    assert ocp._compiled.solver is s2
