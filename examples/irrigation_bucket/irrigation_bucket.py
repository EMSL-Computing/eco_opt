"""Root-zone irrigation bucket, inspired by Shang et al. (2020) water balance.

Not a reproduction of their robust MPC or weather product. Continuous-time
bucket, dry forecast in the OCP, sinusoidal ET, optional rain only on the
plant::

    dS/dt = u + P(t) - E(t) - k S
    min ∫ u dt
    S_min <= S <= S_max,  0 <= u <= u_max,  S(0)=S0,  T fixed.

S is root-zone storage (mm), u irrigation (mm d^{-1}).
"""

import casadi as cs

import eco_opt
from eco_opt import Namespace

T_HORIZON = 14.0
S0 = 40.0
S_MIN = 20.0
S_MAX = 80.0
K_DRAIN = 0.05
E0 = 4.0
E_AMP = 0.6
U_MAX = 15.0
P_RAIN = 10.0
PULSE_T0 = 5.0
PULSE_TF = 7.0

N_INTERVALS = 28
DEGREE = 3

# Frozen after the default-mesh dry open-loop solve (N=28, degree 3, Radau).
REGRESSION_OBJECTIVE = 51.72885336135217


def et(t):
    """Daily sinusoidal ET (mm d^{-1}); always positive."""
    return E0 * (1.0 + E_AMP * cs.sin(2.0 * cs.pi * t))


def rain_pulse(t):
    return P_RAIN * cs.if_else(
        cs.logic_and(t >= PULSE_T0, t < PULSE_TF), 1.0, 0.0
    )


def dry_model(t, x, u):
    S, irr = x[0], u[0]
    E = et(t)
    dS = irr - E - K_DRAIN * S
    y = Namespace(S=S, u=irr, E=E, P=0.0)
    return dS, y


def rain_plant(t, x, u):
    S, irr = x[0], u[0]
    E = et(t)
    P = rain_pulse(t)
    dS = irr + P - E - K_DRAIN * S
    y = Namespace(S=S, u=irr, E=E, P=P)
    return dS, y


def _constraints(sys):
    return [
        (S_MIN, "<=", sys.x[0], "<=", S_MAX),
        (0.0, "<=", sys.u[0], "<=", U_MAX),
        (0.0, "==", eco_opt.t_0(sys.t)),
        (T_HORIZON, "==", eco_opt.t_f(sys.t)),
        (S0, "==", eco_opt.t_0(sys.x[0])),
    ]


def build_ocp():
    sys = eco_opt.System(states=1, controls=1, model=dry_model, name="bucket")
    ocp = eco_opt.Ocp()
    ocp.min([sys.y.u])
    ocp.st(system=sys, constraints=_constraints(sys))
    return ocp, sys


def build_rain_plant():
    return eco_opt.System(states=1, controls=1, model=rain_plant, name="bucket_rain")


def solve(control_intervals=N_INTERVALS, polynomial_degree=DEGREE):
    ocp, sys = build_ocp()
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="radau",
        initial_guess={sys.x[0]: S0, sys.u[0]: E0},
        tf_guess=T_HORIZON,
        ipopt_options={"max_iter": 3000, "print_level": 0, "sb": "yes"},
    )
    return sol, sys, ocp


if __name__ == "__main__":
    sol, sys, _ = solve()
    print(sol)
    print("water use ∫u dt:", sol.objective)
    print("S(T):", sol.x[-1, 0], "min S:", sol.x[:, 0].min())
