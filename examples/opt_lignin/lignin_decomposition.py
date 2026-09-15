"""Optimal ligninolytic activity during litter decomposition
(Chakrawal et al. 2024, OptLignin, https://doi.org/10.5281/zenodo.10553503).

Problem (fixed terminal time)::

    maximize   int_0^TT G dt
    subject to CTdot = -Dh - Do + G,   CT(0) = CT_0
               Codot = -Do + mo*G,     Co(0) = Co_0
               0 <= CT,  0 <= Co
               0 <= vo <= vomax
               vh_max, mo, ro fixed parameters (box-equality "parameters")
               t_f == TT (fixed final time)

where:
    L  = Co / CT                          (aromatic/lignin fraction remaining)
    Dh = vh_max * g(L) * (CT - Co)        (non-aromatic ["holocellulose"] decomposition)
    Do = vo * Co                          (aromatic/lignin-targeting decomposition, control)
    e  = emax - ro * vo                   (growth efficiency, decreasing in effort vo)
    G  = e * (Dh + Do)                    (growth/objective integrand)
    g(L) = exp(-(L/a)^b)                  (decomposability-vs-aromaticity kernel)

``CT`` is total remaining litter carbon, ``Co`` the aromatic (lignin-derived)
carbon fraction of it, ``vo`` the control (rate of investment in
lignin-degrading -- "aromatic-targeting" -- enzymes), maximizing cumulative
microbial growth ``G`` over a fixed decomposition horizon ``TT``.

Fixed kinetic parameters are symbolic ``sys.p`` components pinned by
box-equality constraints. Each solve builds a fresh NLP.

Reference
---------
No single closed-form objective value is published; the paper's Figure 2
instead shows *qualitative* sensitivity of the CT/Co/vo trajectories to
sweeps over (vh_max, ro, ARC0, mo). We reproduce that same sweep in
``plot_figure2.py`` and validate quantitatively that (a) the direct
collocation converges (`Solve_Succeeded`) for every parameter combination,
(b) mass-balance sanity holds (``CT`` and ``Co`` remain non-negative,
``Co <= CT``), and (c) a baseline configuration's objective is
regression-tracked to catch silent transcription regressions.
"""

import casadi as cs

import eco_opt
from eco_opt import Namespace

# g-function parameters (Chakrawal et al. 2024, "exp" scenario)
A = 0.1537
B = 2.0


def g_kernel(L, a=A, b=B):
    return cs.exp(-((L / a) ** b))


def emax_fun(CN0):
    """Maximum growth efficiency as a function of litter C:N ratio."""
    return cs.fmin(6.25 * CN0 ** -0.77, 0.4)


def lignin_model(t, x, u, p):
    CT, Co = x[0], x[1]
    vo = u[0]
    vh_max, mo, ro = p[0], p[1], p[2]

    L = Co / CT
    Dh = vh_max * g_kernel(L) * (CT - Co)
    Do = vo * Co
    e = EMAX - ro * vo
    G = e * (Dh + Do)

    dCT = -Dh - Do + G
    dCo = -Do + mo * G

    y = Namespace(CT=CT, Co=Co, vo=vo, Dh=Dh, Do=Do, G=G, L=L)
    return [dCT, dCo], y


# Baseline parameters (CN0=50)
CN0_BASELINE = 50.0
EMAX = float(6.25 * CN0_BASELINE ** -0.77 if 6.25 * CN0_BASELINE ** -0.77 < 0.4 else 0.4)
VOMAX = 0.1
F = 1.0
MO_DEFAULT = 0.1

VH_LOOP = 0.0015
RO_LOOP = 25.0
ARC0_LOOP = 0.05
CT_0_DEFAULT = 0.5
TT_DEFAULT = 2000.0


def build_ocp(vh_max, mo, ro, CT_0=CT_0_DEFAULT, ARC0=ARC0_LOOP, TT=TT_DEFAULT):
    Co_0 = CT_0 * ARC0

    sys = eco_opt.System(states=2, controls=1, parameters=3, model=lignin_model, name="lignin")
    time = sys.t
    lignin = sys.y

    ocp = eco_opt.Ocp()
    ocp.max([lignin.G])  # Lagrange: maximize integral of G
    ocp.st(
        system=sys,
        constraints=[
            (0, "<=", lignin.CT),
            (0, "<=", lignin.Co),
            (0, "<=", lignin.vo, "<=", VOMAX),
            (vh_max, "==", sys.p[0]),
            (mo, "==", sys.p[1]),
            (ro, "==", sys.p[2]),
            (0, "==", eco_opt.t_0(time)),
            (CT_0, "==", eco_opt.t_0(lignin.CT)),
            (Co_0, "==", eco_opt.t_0(lignin.Co)),
            (TT, "==", eco_opt.t_f(time)),
        ],
    )
    return ocp, sys


def apply_scenario(ocp, sys, vh_max, mo, ro, CT_0=CT_0_DEFAULT, ARC0=ARC0_LOOP):
    """Retune pinned kinetics and initial pools on a compiled NLP."""
    ocp.set_constraints(
        [
            (vh_max, "==", sys.p[0]),
            (mo, "==", sys.p[1]),
            (ro, "==", sys.p[2]),
            (CT_0, "==", eco_opt.t_0(sys.y.CT)),
            (CT_0 * ARC0, "==", eco_opt.t_0(sys.y.Co)),
        ]
    )
    return ocp


def solve(vh_max=VH_LOOP, mo=MO_DEFAULT, ro=RO_LOOP, CT_0=CT_0_DEFAULT, ARC0=ARC0_LOOP,
          TT=TT_DEFAULT, control_intervals=30, polynomial_degree=3):
    ocp, sys = build_ocp(vh_max, mo, ro, CT_0=CT_0, ARC0=ARC0, TT=TT)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        ipopt_options={"max_iter": 2000},
        tf_guess=TT,
    )
    return sol, sys


def compile_ocp(vh_max=VH_LOOP, mo=MO_DEFAULT, ro=RO_LOOP, CT_0=CT_0_DEFAULT,
                ARC0=ARC0_LOOP, TT=TT_DEFAULT, control_intervals=30,
                polynomial_degree=3):
    """Build the lignin NLP once for a parameter sweep."""
    ocp, sys = build_ocp(vh_max, mo, ro, CT_0=CT_0, ARC0=ARC0, TT=TT)
    ocp.build(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        tf_guess=TT,
    )
    return ocp, sys


# Regression-tracked value: no external closed-form objective available (see
# module docstring); computed with eco_opt at the baseline parameters above,
# control_intervals=30, polynomial_degree=3, legendre collocation.
REGRESSION_OBJECTIVE = 0.16281350935527988


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective (integral of G):", sol.objective)
    print("CT(tf)/CT(0):", sol.x[-1, 0] / sol.x[0, 0])
    print("Co(tf)/Co(0):", sol.x[-1, 1] / sol.x[0, 1])
