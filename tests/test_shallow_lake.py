"""Shallow lake: hysteresis of f, two management locals, restoration."""
import os
import sys as _sys

import numpy as np

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "shallow_lake"))

import shallow_lake as lake


def test_three_equilibria_in_hysteresis_window():
    eq = lake.equilibria_at_loading(0.05)
    assert eq.size == 3
    assert eq[0] < 0.3 < eq[1] < 0.9 < eq[2]
    # Middle branch is unstable: dP/dt changes sign through it at this L.
    L = 0.05
    def dP(P):
        return L - lake.B * P + lake.recycling(P)

    assert dP(eq[0] + 0.02) < 0  # toward clear
    assert dP(eq[2] - 0.02) > 0  # toward turbid


def test_management_guess_splits_two_locals():
    sol_lo, _ = lake.solve_management(L_guess=0.05)
    sol_hi, _ = lake.solve_management(L_guess=0.50)
    assert sol_lo.success and sol_hi.success
    assert abs(sol_lo.objective - lake.REGRESSION_J_LOW_GUESS) < 0.05
    assert abs(sol_hi.objective - lake.REGRESSION_J_HIGH_GUESS) < 0.05
    # High loading guess finds the better (less negative) welfare and a
    # cleaner terminal lake — two locals, not noise.
    assert sol_hi.objective > sol_lo.objective + 1.0
    assert sol_hi.x[-1, 0] < sol_lo.x[-1, 0] - 0.3


def test_restoration_reaches_clear_target():
    sol, _ = lake.solve_restoration()
    assert sol.success, sol.return_status
    assert sol.x[-1, 0] <= lake.P_CLEAN + 1e-4
    assert sol.x[0, 0] == lake.P0_RESTORE
    # Must cut below the desired long-run loading to leave the turbid well.
    assert sol.u[:, 0].min() < lake.L_DESIRED - 1e-3
