"""Analytic-reference tests for Cohen (1971) and King & Roughgarden (1982)."""
import os
import sys as _sys

import numpy as np

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "cohen_king_roughgarden"))

import cohen_terminal_yield as cohen
import king_singular_allocation as king


def test_cohen_matches_analytic_bang_bang():
    sol, sys = cohen.solve()
    assert sol.success, sol.return_status
    assert abs(sol.objective - cohen.ANALYTIC_YIELD) < 5e-2
    tF = cohen.flowering_time_from_ocp(sol)
    assert abs(tF - cohen.ANALYTIC_SWITCH) < 0.15
    # bang-bang: vegetative then reproductive
    assert sol.u[0, 0] > 0.85
    assert sol.u[-1, 0] < 0.15
    # Cohen corollary: equal conversion efficiencies ⇒ x1(T) = x2(T)
    assert abs(sol.x[-1, 0] - sol.x[-1, 1]) < 0.5
    assert abs(sol.x[-1, 1] - cohen.ANALYTIC_YIELD) < 5e-2
    assert abs(sol.x[0, 0] - cohen.X10) < 1e-8
    assert abs(sol.x[0, 1]) < 1e-6


def test_king_junction_constants_match_published_values():
    w, y, z = king.junction_constants()
    assert abs(y - 2.79328213) < 1e-8
    assert abs(z - 0.55763674) < 1e-8
    assert abs(y + z - 3.35091887) < 1e-8
    # residual of e^w = w^2 + w + 1
    assert abs(np.exp(w) - (w * w + w + 1.0)) < 1e-12


def test_king_matches_analytic_singular_manifold():
    sol, sys = king.solve()
    assert sol.success, sol.return_status
    assert abs(sol.objective - king.ANALYTIC_OBJECTIVE) < 0.15
    # Terminal reproductive/vegetative ratio is the robust junction diagnostic.
    term = sol.x[-1, 1] / sol.x[-1, 0]
    assert abs(term / king.ANALYTIC_TERMINAL_RATIO - 1.0) < 0.08
    t2 = king.t2_from_control(sol)
    assert abs(t2 - king.ANALYTIC_T2) < 0.15
    # On the interior of the singular arc the state-ratio test holds even
    # if collocation chatters in u(t).
    t = sol.time
    interior = (t > 0.25) & (t < king.ANALYTIC_T2 - 0.25)
    ratio = king.singular_ratio(t[interior], sol.x[interior, 0], sol.x[interior, 1])
    assert np.mean(np.abs(ratio - 1.0)) < 0.12
    assert sol.x[0, 0] == king.X10 or abs(sol.x[0, 0] - king.X10) < 1e-8
    assert abs(sol.x[0, 1] - king.X20) < 1e-8
    # After t2 the control is the reproductive bang.
    t_u = sol.time_intervals[:-1]
    assert np.mean(sol.u[t_u > king.ANALYTIC_T2 + 0.15, 0]) < 0.2
