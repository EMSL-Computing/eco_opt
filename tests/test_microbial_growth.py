"""Regression + physical-sanity test for the microbial growth OCP
(Manzoni et al. 2022, see examples/microbial_growth/README.md).
"""
import os
import sys as _sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE_DIR = os.path.join(HERE, "examples", "microbial_growth")


def _load_example():
    if EXAMPLE_DIR not in _sys.path:
        _sys.path.insert(0, EXAMPLE_DIR)
    import microbial_growth as mg

    return mg


def test_microbial_growth_matches_regression_and_is_physical():
    mg = _load_example()
    sol, sys = mg.solve()

    assert sol.success, f"IPOPT did not report success: {sol.return_status}"

    # -- physical sense: free terminal time driven by resource exhaustion --
    assert sol.x[-1, 0] < 1e-4, f"C(T)={sol.x[-1,0]} should be ~0 (resource exhausted)"
    assert abs(sol.x[0, 0] - 1.0) < 1e-6, "C(0) should equal 1"

    # -- box feasibility ----------------------------------------------------
    assert np.all(sol.x[:, 0] >= mg.C_MIN - 1e-6)
    assert np.all(sol.x[:, 0] <= mg.C_MAX + 1e-6)
    assert np.all(sol.u[:, 0] >= mg.U_MIN - 1e-6)
    assert np.all(sol.u[:, 0] <= mg.U_MAX + 1e-6)

    # -- regression tracking --------------------------------------------------
    assert abs(sol.objective - mg.REGRESSION_OBJECTIVE) < 1e-4
    assert abs(sol.tf - mg.REGRESSION_TF) < 1e-3

    t_an, C_an, u_an, _g_an, T_an, J_an = mg.analytical_solution()
    assert abs(T_an - sol.tf) < 5e-3
    assert abs(J_an - sol.objective) < 5e-4
    u_col = np.interp(t_an, sol.time[:-1], sol.u[:, 0])
    scale = np.max(np.abs(u_an))
    assert np.max(np.abs(u_col - u_an)) < 0.05 * scale
