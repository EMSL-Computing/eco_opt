"""Analytic-reference tests for Iwasa–Roughgarden Liebig shoot/root/fruit."""
import os
import sys as _sys

import numpy as np

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "iwasa_roughgarden_shoot_root"))

import shoot_root as ir


def test_iwasa_liebig_closed_forms():
    assert abs(ir.GAMMA - 0.5) < 1e-12
    assert abs(ir.ANALYTIC_T_S - 4.0) < 1e-12
    assert abs(ir.ANALYTIC_T_BAL - np.log(5.0)) < 1e-12
    expected_fruit = 2.0 * np.exp(2.0) / np.sqrt(5.0)
    assert abs(ir.ANALYTIC_FRUIT - expected_fruit) < 1e-10
    u_S, u_R = ir.balanced_controls()
    assert abs(u_S - 0.5) < 1e-12
    assert abs(u_R - 0.5) < 1e-12


def test_iwasa_ocp_recovers_three_phases():
    sol, sys = ir.solve()
    assert sol.success, sol.return_status
    # Terminal fruit against the true Liebig closed form (smooth-min transcription).
    assert abs(sol.x[-1, 2] / ir.ANALYTIC_FRUIT - 1.0) < 0.03
    assert abs(sol.objective / ir.ANALYTIC_FRUIT - 1.0) < 0.03

    t_bal = ir.balance_time_from_ocp(sol)
    t_s = ir.fruit_switch_from_ocp(sol)
    assert abs(t_bal - ir.ANALYTIC_T_BAL) < 0.25
    assert abs(t_s - ir.ANALYTIC_T_S) < 0.25
    assert t_bal < t_s < ir.T_SEASON

    # Phase 1 (root-deficient): almost all photosynthate to root.
    t_u = sol.time_intervals[:-1]
    early = t_u < ir.ANALYTIC_T_BAL - 0.2
    if np.any(early):
        assert np.mean(sol.u[early, 1]) > 0.7
        assert np.mean(sol.u[early, 0]) < 0.3

    # Balanced window: R/S near α/β = 1, both organs growing.
    t = sol.time
    balanced = (t > ir.ANALYTIC_T_BAL + 0.25) & (t < ir.ANALYTIC_T_S - 0.25)
    ratio = sol.x[balanced, 1] / np.maximum(sol.x[balanced, 0], 1e-12)
    assert np.mean(np.abs(ratio - ir.ALPHA / ir.BETA)) < 0.2

    # After fruit switch, S and R freeze (no further vegetative growth).
    late = t >= ir.ANALYTIC_T_S + 0.2
    if np.sum(late) >= 2:
        assert np.ptp(sol.x[late, 0]) < 0.15 * sol.x[late, 0][0]
        assert np.ptp(sol.x[late, 1]) < 0.15 * sol.x[late, 1][0]
        # Fruit keeps accumulating.
        assert sol.x[-1, 2] > sol.x[late, 2][0] + 0.5

    # Simplex path constraint.
    veg = sol.u[:, 0] + sol.u[:, 1]
    assert np.all(veg >= -1e-4)
    assert np.all(veg <= 1.0 + 1e-3)
    assert abs(sol.x[0, 0] - ir.S0) < 1e-8
    assert abs(sol.x[0, 1] - ir.R0) < 1e-8
    assert abs(sol.x[0, 2] - ir.F0) < 1e-8
