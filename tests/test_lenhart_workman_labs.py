"""Lenhart and Workman labs: analytic facts from the book."""
import os
import sys as _sys

import numpy as np

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "lenhart_workman_labs"))

import lab1_introductory as lab1
import lab2_mold_fungicide as lab2
import lab3_bacteria as lab3
import lab4_bounded as lab4
import lab5_cancer as lab5
import lab6_fish as lab6
import lab7_epidemic as lab7
import lab8_hiv as lab8
import lab9_bears as lab9
import lab10_glucose as lab10
import lab11_timber as lab11
import lab12_bioreactor as lab12
import lab13_predator_prey as lab13


def test_lab1_default_succeeds_and_front_loads_control():
    sol, _sys = lab1.solve()
    assert sol.success
    assert abs(sol.tf - lab1.T_HORIZON) < 1e-8
    assert abs(sol.x[0, 0] - lab1.DEFAULT["x0"]) < 1e-6
    # Book Fig. 5.1: u starts high and falls toward 0.
    assert sol.u[0, 0] > sol.u[-1, 0]
    assert sol.u[-1, 0] < 0.3


def test_lab1_A_zero_implies_zero_control():
    sol, _sys = lab1.solve(A=0.0)
    assert sol.success
    assert np.max(np.abs(sol.u[:, 0])) < 1e-4


def test_lab1_trajectories_depend_on_A_over_B_only():
    sol_a, _ = lab1.solve(A=1.0, B=2.0)
    sol_b, _ = lab1.solve(A=2.0, B=4.0)
    assert sol_a.success and sol_b.success
    t = sol_a.time
    xb = np.interp(t, sol_b.time, sol_b.x[:, 0])
    ub = np.interp(sol_a.time_intervals[:-1], sol_b.time_intervals[:-1], sol_b.u[:, 0])
    assert np.sqrt(np.mean((sol_a.x[:, 0] - xb) ** 2)) < 1e-3
    assert np.sqrt(np.mean((sol_a.u[:, 0] - ub) ** 2)) < 1e-3


def test_lab2_default_control_ends_near_zero():
    sol, _sys, p = lab2.solve()
    assert sol.success
    assert abs(sol.tf - p["T"]) < 1e-8
    assert abs(sol.x[0, 0] - p["x0"]) < 1e-6
    assert np.all(sol.u[:, 0] > -1e-4)
    # Book: u rises, plateaus, then declines toward the end.
    mid = sol.u[len(sol.u) // 2, 0]
    assert mid > sol.u[0, 0]
    assert sol.u[-1, 0] < 0.5 * mid


def test_lab2_A10_holds_mold_below_untreated():
    sol, _sys, p = lab2.solve(**lab2.CASE_A10)
    assert sol.success
    x_free = lab2.untreated_state(sol.time, p["r"], p["M"], p["x0"])
    # Book Fig. 6.1: treated mold is held below the u≡0 trajectory.
    assert np.all(sol.x[:, 0] < x_free + 1e-6)
    n = len(sol.x)
    mid = sol.x[n // 2, 0]
    assert mid < p["x0"]
    assert sol.x[-1, 0] > mid
    umid = sol.u[len(sol.u) // 2, 0]
    assert sol.u[-1, 0] < 0.5 * umid


def test_lab3_default_back_loads_nutrient():
    sol, _sys, p = lab3.solve()
    assert sol.success
    assert abs(sol.x[0, 0] - p["x0"]) < 1e-6
    # Book: nutrient is concentrated at the end of [0, 1].
    assert sol.u[-1, 0] > sol.u[0, 0]
    assert sol.x[-1, 0] > sol.x[0, 0]


def test_lab3_A_zero_implies_zero_control():
    sol, _sys, _p = lab3.solve(A=0.0)
    assert sol.success
    assert np.max(np.abs(sol.u[:, 0])) < 1e-3


def test_lab4_wide_bounds_match_lab1():
    sol4, _, p = lab4.solve()
    sol1, _ = lab1.solve(A=p["A"], B=1.0, C=p["C"], x0=p["x0"])
    assert sol4.success
    u1 = np.interp(sol4.time_intervals[:-1], sol1.time_intervals[:-1], sol1.u[:, 0])
    assert np.sqrt(np.mean((sol4.u[:, 0] - u1) ** 2)) < 5e-2
    assert np.all(sol4.u[:, 0] >= p["M1"] - 1e-8)
    assert np.all(sol4.u[:, 0] <= p["M2"] + 1e-8)


def test_lab4_equal_bounds_force_that_value():
    sol, _, _ = lab4.solve(M1=0.0, M2=0.0)
    assert sol.success
    assert np.max(np.abs(sol.u[:, 0])) < 1e-6


def test_lab5_treated_below_untreated_gompertz():
    sol, _, p = lab5.solve()
    assert sol.success
    assert np.all(sol.u[:, 0] >= -1e-8)
    Nfree = lab5.untreated_N(sol.time, p["r"], p["N0"])
    assert np.all(sol.x[:, 0] < Nfree + 1e-4)
    assert sol.u[0, 0] > sol.u[-1, 0]


def test_lab6_harvest_starts_late_and_Ak_product():
    sol, _, p = lab6.solve()
    assert sol.success
    assert sol.u[0, 0] < 0.15
    assert sol.u.max() > sol.u[0, 0] + 0.1
    sol_a, _, _ = lab6.solve(A=10.0, k=5.0)
    sol_b, _, _ = lab6.solve(A=5.0, k=10.0)
    ua = np.interp(sol_b.time_intervals[:-1], sol_a.time_intervals[:-1], sol_a.u[:, 0])
    assert np.sqrt(np.mean((sol_b.u[:, 0] - ua) ** 2)) < 5e-2


def test_lab7_vaccination_in_bounds_and_reduces_I():
    sol, _, p = lab7.solve()
    assert sol.success
    assert np.all(sol.u[:, 0] >= -1e-8)
    assert np.all(sol.u[:, 0] <= 0.9 + 1e-8)
    assert sol.x[-1, 2] < sol.x[0, 2]


def test_lab8_therapy_strongest_first():
    sol, _, _ = lab8.solve()
    assert sol.success
    # u=0 is maximum therapy; book starts at strongest dose.
    assert sol.u[0, 0] < sol.u[-1, 0]
    assert sol.x[-1, 0] > sol.x[0, 0]


def test_lab9_park_harvest_near_zero():
    sol, _, _ = lab9.solve()
    assert sol.success
    assert np.max(sol.u[:, 0]) < 0.05
    assert np.max(sol.u[:, 1]) > np.max(sol.u[:, 0])


def test_lab10_glucose_moves_toward_target():
    sol, _, p = lab10.solve()
    assert sol.success
    assert abs(sol.x[-1, 0] - p["l"]) < abs(sol.x[0, 0] - p["l"])


def test_lab11_bang_bang_switch_matches_analytic():
    sol, _, p = lab11.solve()
    assert sol.success
    t_star = lab11.analytic_switch(p["k"], p["r"], p["T"])
    t = sol.time_intervals[:-1]
    u = sol.u[:, 0]
    idx = np.argmax(u < 0.5)
    assert abs(t[idx] - t_star) < 0.15
    assert u[0] > 0.8
    assert u[-1] < 0.2


def test_lab12_K2_no_injection_K21_injects():
    sol0, _, _ = lab12.solve()
    assert sol0.success
    assert np.mean(sol0.u[:, 0]) < 0.15
    sol1, _, p = lab12.solve(**lab12.INJECT)
    assert sol1.success
    assert sol1.u[0, 0] > 0.8 * p["M"]


def test_lab13_isoperimetric_budget():
    sol, _, p = lab13.solve()
    assert sol.success
    assert abs(sol.x[-1, 2] - p["B"]) < 1e-3
    assert np.all(sol.u[:, 0] >= -1e-8)
    assert np.all(sol.u[:, 0] <= p["M"] + 1e-8)
