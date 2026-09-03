"""Irrigation bucket: dry open-loop, same-plant closed loop, rain plant."""
import os
import sys as _sys

import numpy as np

import eco_opt

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "irrigation_bucket"))

import irrigation_bucket as ib
import open_vs_closed as ovc


def test_dry_open_loop_hits_floor_and_regression_J():
    sol, sys, _ = ib.solve()
    assert sol.success, sol.return_status
    assert abs(sol.tf - ib.T_HORIZON) < 1e-10
    assert abs(sol.x[0, 0] - ib.S0) < 1e-8
    assert sol.x[:, 0].min() >= ib.S_MIN - 1e-4
    assert sol.x[:, 0].max() <= ib.S_MAX + 1e-4
    assert np.all(sol.u[:, 0] >= -1e-8)
    assert np.all(sol.u[:, 0] <= ib.U_MAX + 1e-8)
    # Rides the floor: without irrigation S would drain below S_min.
    assert sol.x[-1, 0] < ib.S_MIN + 0.05
    assert abs(sol.objective - ib.REGRESSION_OBJECTIVE) < 1e-4


def test_shrinking_closed_loop_matches_dry_open_loop():
    sol, sys, _ = ib.solve()
    cl = ovc.closed_loop(plant=None)
    assert cl.success
    assert np.max(np.abs(cl.u[:, 0] - sol.u[:, 0])) < 0.08
    assert np.max(np.abs(cl.x[:, 0] - sol.x[:, 0])) < 0.5


def test_closed_loop_irrigates_less_after_unforecast_rain():
    sol, sys, ocp = ib.solve()
    plant = ib.build_rain_plant()
    ol_rain = eco_opt.simulate(plant, [ib.S0], sol.time, u=sol)
    cl_rain = ovc.closed_loop(plant=plant)
    assert cl_rain.success
    J_tape = ovc.water_use(ol_rain.time, ol_rain.u[:, 0])
    J_cl = ovc.water_use(cl_rain.time, cl_rain.u[:, 0])
    # Same u* replayed: same water. Closed loop cuts irrigation after the pulse.
    assert abs(J_tape - sol.objective) < 0.2
    assert J_cl < J_tape - 5.0
    # Tape keeps watering a wetter soil; closed loop does not pile onto S_max.
    assert cl_rain.x[:, 0].max() <= ol_rain.x[:, 0].max() + 1e-6
