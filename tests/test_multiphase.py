"""Multi-phase transcription: Cohen two-phase bang-bang with free switch."""
import os
import sys as _sys

import numpy as np

import eco_opt

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "multiphase_cohen"))

import cohen_two_phase as cohen


def test_phase_is_public():
    assert hasattr(eco_opt, "Phase")
    p = eco_opt.Phase(eco_opt.System(states=1, controls=0, model=lambda t, x, u: x), constraints=[])
    assert p.system.nx == 1


def test_cohen_two_phase_matches_analytic():
    sol, systems = cohen.solve()
    assert sol.success, sol.return_status
    assert len(sol.phases) == 2
    assert abs(sol.t0) < 1e-12
    assert abs(sol.tf - cohen.T_SEASON) < 1e-8
    assert abs(sol.switch_times[0] - cohen.ANALYTIC_SWITCH) < 1e-5
    assert abs(sol.objective - cohen.ANALYTIC_YIELD) < 1e-6
    # junction continuity
    assert abs(sol.phases[0].tf - sol.phases[1].t0) < 1e-8
    assert np.allclose(sol.phases[0].x[-1], sol.phases[1].x[0], atol=1e-7)
    # Cohen corollary: x1(T) = x2(T) = exp(T-1)
    assert abs(sol.x[-1, 0] - cohen.ANALYTIC_YIELD) < 1e-5
    assert abs(sol.x[-1, 1] - cohen.ANALYTIC_YIELD) < 1e-5
    # vegetative phase: repro stays ~0; reproductive phase: veg frozen
    assert np.max(np.abs(sol.phases[0].x[:, 1])) < 1e-6
    assert np.ptp(sol.phases[1].x[:, 0]) < 1e-5


def test_cohen_two_phase_recovers_switch_from_bad_knot():
    sol, systems = cohen.solve(tf_guess=[0.0, 2.0, 5.0])
    assert sol.success
    assert abs(sol.switch_times[0] - cohen.ANALYTIC_SWITCH) < 1e-5
    assert abs(sol.objective - cohen.ANALYTIC_YIELD) < 1e-6


def test_cohen_two_phase_optimize_matches_solve():
    sol_s, _ = cohen.solve()
    sol_o, _ = cohen.optimize()
    assert sol_s.success and sol_o.success
    assert abs(sol_s.objective - sol_o.objective) < 1e-8
    assert abs(sol_s.switch_times[0] - sol_o.switch_times[0]) < 1e-8
