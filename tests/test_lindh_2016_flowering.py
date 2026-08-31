"""Regression/validation tests for the Lindh et al. (2016) flowering-time
reproduction (see examples/lindh_2016_flowering/README.md).

Validates the paper's two central claims:
  1. (Theorem 1) the optimal control is bang-bang for all three growth
     models (exponential, logistic, WBE).
  2. (Fig. 2 g-i) the direction of d(tF*)/dP flips with growth
     constraint: monotone increasing for exponential growth, but
     hump-shaped (interior maximum, i.e. NOT monotone) for logistic and
     WBE growth.
"""
import os
import sys as _sys

import numpy as np

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "lindh_2016_flowering"))

from optimal_flowering import (
    solve_ocp, is_bang_bang, flowering_time_from_ocp,
    optimal_flowering_time_1d,
)


def test_bang_bang_control_for_all_growth_models():
    for growth in ["exponential", "logistic", "wbe"]:
        sol, sys = solve_ocp(growth, P=2.0)
        assert sol.success, f"{growth} OCP failed to converge"
        assert is_bang_bang(sol), f"{growth}: optimal control is not bang-bang"


def test_ocp_and_1d_flowering_time_agree():
    """Cross-validate the full-OCP switch time against the fast 1-D
    closed-form/scalar-optimization approach (independent methods)."""
    for growth in ["exponential", "logistic", "wbe"]:
        sol, sys = solve_ocp(growth, P=2.0)
        tF_ocp = flowering_time_from_ocp(sol)
        tF_1d = optimal_flowering_time_1d(growth, P=2.0)
        assert abs(tF_ocp - tF_1d) < 0.05, (
            f"{growth}: OCP tF={tF_ocp:.4f} vs 1D tF={tF_1d:.4f} disagree"
        )


def test_exponential_growth_monotone_increasing_flowering_time():
    P_sweep = np.linspace(0.3, 15.0, 40)
    tF = np.array([optimal_flowering_time_1d("exponential", P) for P in P_sweep])
    assert np.all(np.diff(tF) >= -1e-6), "exponential tF*(P) should be monotone increasing"


def test_constrained_growth_flips_direction_interior_maximum():
    """The paper's central claim: constrained growth (logistic, WBE) has
    an INTERIOR MAXIMUM in tF*(P), i.e. tF* eventually DECREASES as P
    keeps increasing -- unlike unconstrained (exponential) growth."""
    P_sweep = np.linspace(0.3, 15.0, 40)
    for growth in ["logistic", "wbe"]:
        tF = np.array([optimal_flowering_time_1d(growth, P) for P in P_sweep])
        argmax_idx = int(np.argmax(tF))
        assert 0 < argmax_idx < len(tF) - 1, (
            f"{growth}: expected an interior maximum in tF*(P), got argmax at "
            f"boundary index {argmax_idx}/{len(tF)-1}"
        )
        # after the peak, tF* must decrease (flip in direction)
        assert tF[-1] < tF[argmax_idx] - 1e-3, (
            f"{growth}: tF*(P) did not decrease after its interior maximum"
        )
