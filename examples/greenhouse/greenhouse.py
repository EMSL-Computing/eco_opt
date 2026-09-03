"""Greenhouse climate control with tabulated light and outside temperature.

van Straten, van Ooteghem, van Willigenburg & van Henten, *Optimal
Control of Greenhouse Cultivation*. Crop dry weight x1 and indoor
temperature x2; heating u is the control. Light I(t) and outside
temperature T_out(t) are exogenous (not controls): two-day sinusoids,
clipped so I >= 0, stored as linear interpolants of t.

    minimize   -p5 * x1(T) + int_0^T p4 * u dt
    subject to dx1/dt = p1 * I(t) * x2
               dx2/dt = p2 * (T_out(t) - x2) + p3 * u
               x1(0)=0, x2(0)=10, T=48 h fixed
               0 <= u <= 10

Profit is crop value minus heating cost (min of the negative profit).
No closed form; J is regression-tracked.
"""

import numpy as np
import casadi as cs

import eco_opt
from eco_opt import Namespace

T_HORIZON = 48.0
X1_0 = 0.0
X2_0 = 10.0
U_MIN, U_MAX = 0.0, 10.0

P1 = 7.5e-8
P2 = 1.0
P3 = 0.1
P4 = 4.55e-4
P5 = 136.4

# Table used by the in-code interpolants (same sampling as the published listing).
_TE = np.arange(-1.0, 49.0 + 1e-12, 0.2)
_I_TAB = np.maximum(0.0, 800.0 * np.sin(4.0 * np.pi * _TE / T_HORIZON - 0.65 * np.pi))
_TOUT_TAB = 15.0 + 10.0 * np.sin(4.0 * np.pi * _TE / T_HORIZON - 0.65 * np.pi)


def interpolants_from_arrays(t_grid, I_tab, T_out_tab, tag="inline"):
    """Linear interpolants of light and outside temperature vs time."""
    t_grid = np.asarray(t_grid, dtype=float).ravel()
    I_tab = np.asarray(I_tab, dtype=float).ravel()
    T_out_tab = np.asarray(T_out_tab, dtype=float).ravel()
    I_fun = cs.interpolant(f"I_light_{tag}", "linear", [t_grid], I_tab)
    T_fun = cs.interpolant(f"T_out_{tag}", "linear", [t_grid], T_out_tab)
    return I_fun, T_fun


def interpolants_from_file(path, tag="file"):
    """Read CSV with header ``t,light,T_out`` and return interpolants."""
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return interpolants_from_arrays(data[:, 0], data[:, 1], data[:, 2], tag=tag)


def make_model(I_fun, T_fun):
    def greenhouse_model(t, x, u):
        crop, temp = x[0], x[1]
        heat = u[0]
        I = I_fun(t)
        T_out = T_fun(t)
        dx = cs.vertcat(P1 * I * temp, P2 * (T_out - temp) + P3 * heat)
        y = Namespace(crop=crop, temp=temp, heat=heat, light=I, T_out=T_out)
        return dx, y

    return greenhouse_model


_I_FUN, _TOUT_FUN = interpolants_from_arrays(_TE, _I_TAB, _TOUT_TAB, tag="inline")


def light(t):
    """Incoming light I(t) [W m^{-2}]. Works with floats or CasADi t."""
    return _I_FUN(t)


def outside_temp(t):
    """Outside temperature T_out(t) [°C]."""
    return _TOUT_FUN(t)


def greenhouse_model(t, x, u):
    return make_model(_I_FUN, _TOUT_FUN)(t, x, u)


def build_ocp(model=None, name="greenhouse"):
    if model is None:
        model = greenhouse_model
    sys = eco_opt.System(states=2, controls=1, model=model, name=name)
    ocp = eco_opt.Ocp()
    # min -p5 x1(T) + ∫ p4 u dt  ==  max crop value minus heating cost
    ocp.min([eco_opt.t_f(-P5 * sys.y.crop), P4 * sys.y.heat])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (T_HORIZON, "==", eco_opt.t_f(sys.t)),
            (X1_0, "==", eco_opt.t_0(sys.y.crop)),
            (X2_0, "==", eco_opt.t_0(sys.y.temp)),
            (U_MIN, "<=", sys.y.heat, "<=", U_MAX),
        ],
    )
    return ocp, sys


def solve(control_intervals=100, polynomial_degree=5, model=None, name="greenhouse"):
    ocp, sys = build_ocp(model=model, name=name)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: X1_0, sys.x[1]: X2_0, sys.u[0]: 1.0},
        ipopt_options={"max_iter": 2000, "tol": 1e-8},
        tf_guess=T_HORIZON,
    )
    return sol, sys


# Default mesh: N=100, degree 5, Legendre. Frozen after a verified solve.
REGRESSION_OBJECTIVE = -2.826212087862474


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("J =", sol.objective)
    print("crop(T) =", sol.x[-1, 0], "temp(T) =", sol.x[-1, 1])
    print("u min/max", float(sol.u[:, 0].min()), float(sol.u[:, 0].max()))
