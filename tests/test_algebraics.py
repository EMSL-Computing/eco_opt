"""Algebraic states: reduced ODE vs DAE on the QSS microbe toys."""
import os
import sys as _sys

import numpy as np

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "algebraic_qss"))

import linear_qss as lin
import quadratic_qss as quad


def _coarse_u(sol):
    return np.append(sol.u[:, 0], sol.u[-1, 0])


def test_linear_qss_ode_and_dae_match():
    sol_ode, sys_ode = lin.solve_ode()
    sol_dae, sys_dae = lin.solve_dae()
    assert sol_ode.success, sol_ode.return_status
    assert sol_dae.success, sol_dae.return_status
    assert abs(sol_ode.tf - lin.T_HORIZON) < 1e-10
    assert abs(sol_dae.tf - lin.T_HORIZON) < 1e-10
    assert abs(sol_ode.objective - sol_dae.objective) < 1e-4
    assert np.allclose(sol_ode.x[:, 0], sol_dae.x[:, 0], atol=2e-3)
    assert np.allclose(sol_ode.u[:, 0], sol_dae.u[:, 0], atol=2e-3)

    # DAE M matches the explicit formula on the coarse grid.
    u = _coarse_u(sol_dae)
    M_formula = np.array(
        [float(lin.biomass(c, ui)) for c, ui in zip(sol_dae.x[:, 0], u)]
    )
    assert np.allclose(sol_dae.z[:, 0], M_formula, atol=5e-3)
    assert np.all(sol_dae.z[:, 0] >= -1e-8)
    assert np.all(sol_dae.u[:, 0] >= -1e-8)
    assert np.all(sol_dae.u[:, 0] <= lin.UMAX + 1e-8)


def test_quadratic_qss_dae_matches_positive_root():
    sol_ode, _ = quad.solve_ode()
    sol_dae, _ = quad.solve_dae()
    assert sol_ode.success, sol_ode.return_status
    assert sol_dae.success, sol_dae.return_status
    assert abs(sol_ode.objective - sol_dae.objective) < 2e-3
    assert np.allclose(sol_ode.x[:, 0], sol_dae.x[:, 0], atol=5e-3)

    u = _coarse_u(sol_dae)
    M_formula = np.array(
        [
            float(quad.biomass_positive_root(c, ui))
            for c, ui in zip(sol_dae.x[:, 0], u)
        ]
    )
    assert np.allclose(sol_dae.z[:, 0], M_formula, atol=1e-2)
    # Residual of g at interval starts.
    growth = np.array(
        [float(quad.growth_rate(c, ui)) for c, ui in zip(sol_dae.x[:, 0], u)]
    )
    g = growth - quad.D * sol_dae.z[:, 0] - quad.C_QUAD * sol_dae.z[:, 0] ** 2
    assert np.max(np.abs(g)) < 2e-3
    assert np.all(sol_dae.z[:, 0] >= -1e-8)
