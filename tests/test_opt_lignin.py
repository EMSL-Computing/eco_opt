"""Regression + mass-balance sanity test for the OptLignin Figure 2
litter-decomposition OCP (see examples/opt_lignin/README.md).
"""
import os
import sys as _sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE_DIR = os.path.join(HERE, "examples", "opt_lignin")


def _load_example():
    if EXAMPLE_DIR not in _sys.path:
        _sys.path.insert(0, EXAMPLE_DIR)
    import lignin_decomposition as ld

    return ld


def test_lignin_decomposition_baseline_matches_regression_and_is_physical():
    ld = _load_example()
    sol, sys = ld.solve()

    assert sol.success, f"IPOPT did not report success: {sol.return_status}"

    # -- mass-balance sanity: carbon pools stay non-negative, Co <= CT ------
    assert np.all(sol.x[:, 0] >= -1e-6), "CT went negative"
    assert np.all(sol.x[:, 1] >= -1e-6), "Co went negative"
    assert np.all(sol.x[:, 1] <= sol.x[:, 0] + 1e-6), "Co exceeded CT"

    # -- box feasibility ------------------------------------------------------
    assert np.all(sol.u[:, 0] >= -1e-6)
    assert np.all(sol.u[:, 0] <= ld.VOMAX + 1e-6)

    # -- fixed terminal time honored -----------------------------------------
    assert abs(sol.tf - ld.TT_DEFAULT) < 1e-6

    # -- regression tracking --------------------------------------------------
    assert abs(sol.objective - ld.REGRESSION_OBJECTIVE) < 1e-4


def test_lignin_decomposition_parameter_variant_converges():
    """Sanity: a non-baseline parameter combination (used in the Figure 2
    sweep) also converges and respects mass-balance/box bounds."""
    ld = _load_example()
    sol, sys = ld.solve(vh_max=ld.VH_LOOP * 3, ro=ld.RO_LOOP * 2, mo=0.0,
                         ARC0=ld.ARC0_LOOP * 2)
    assert sol.success, f"IPOPT did not report success: {sol.return_status}"
    assert np.all(sol.x[:, 0] >= -1e-6)
    assert np.all(sol.x[:, 1] >= -1e-6)
