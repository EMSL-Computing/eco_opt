"""King & Roughgarden (1982) graded allocation: bang–singular–bang.

Problem (fixed season, linear vegetative production, log-reproductive
running payoff)::

    maximize   int_0^T ln(x2(t)) dt
    subject to dx1/dt = u * x1
               dx2/dt = (1 - u) * x1
               x1(0) = X10 > 0,  x2(0) = X20 > 0
               0 <= u(t) <= 1,   t in [0, T]

Same two-compartment plant as Cohen (1971); the objective is the
continuous-time geometric-mean / log-seed-yield criterion used by King
and Roughgarden (TPB 22:1–16). Because the running payoff is strictly
increasing and concave in reproductive mass, a mixed-growth *singular*
arc becomes optimal: vegetative and reproductive tissues grow together
for part of the season.

The state equation for ``x2`` is ``(1-u)*x1`` (photosynthate comes from
the vegetative body). Some secondary transcriptions write ``(1-u)*x2``;
that is a typesetting error (Atkins et al., arXiv:2010.06744, HTML eq.
70 vs. the problem statement (72)).

Analytic solution (start on the singular manifold)
--------------------------------------------------
The Hamiltonian is affine in ``u``. Pontryagin therefore yields bang or
singular control. On a singular arc the switching function and its first
two time derivatives vanish, which forces the state-ratio rule

    x1(t) = (T - 1 - t) * x2(t)
    u*(t) = 1 - 1/(T - 1 - t)

The singular arc cannot reach ``t = T`` (the factor ``T-1-t`` would
change sign). It joins a terminal purely-reproductive bang ``u=0`` at a
parameter-free junction. Writing ``y = T - t2`` and
``z = x2(t2)/x1(t2)``, the junction conditions reduce to

    e^w = w^2 + w + 1,   w > 0
    y = 1 + w
    z = 1/w

with numeric values ``y ≈ 2.79328213``, ``z ≈ 0.55763674``, and terminal
ratio ``x2(T)/x1(T) = y + z ≈ 3.35091887`` (Atkins et al. 2020, eqs.
96–99). Hence ``t2 = T - y`` is independent of the initial masses.

If the plant starts *on* the singular manifold,

    X10 = (T - 1) * X20    ⇒    t1 = 0

(the case used here). If ``X10/X20 > T-1``, an initial purely-vegetative
bang ``u=1`` runs until the manifold is hit; that variant is not the
default demo.

Default parameters: T=5, X20=1, X10=4  ⇒  t1=0, t2≈2.2067.

Collocation note: ``u(t)`` on a singular arc is typically polluted by
chatter. Tests and figures therefore treat the *state ratio*
``x1/((T-1-t) x2) ≈ 1`` as the primary diagnostic, not the control.

See ``README.md`` in this directory for the full derivation.
"""

import math

import casadi as cs
import numpy as np

import eco_opt
from eco_opt import Namespace

T_SEASON = 5.0
X20 = 1.0


def junction_w(tol=1e-15, maxiter=80):
    """Positive root of e^w = w^2 + w + 1 (w=0 is the trivial root)."""
    w = 1.8
    for _ in range(maxiter):
        f = math.exp(w) - w * w - w - 1.0
        df = math.exp(w) - 2.0 * w - 1.0
        w_new = w - f / df
        if abs(w_new - w) < tol:
            return w_new
        w = w_new
    return w


def junction_constants():
    """Return ``(w, y, z)`` with y=T-t2, z=x2(t2)/x1(t2)."""
    w = junction_w()
    return w, 1.0 + w, 1.0 / w


def on_manifold_x10(T=T_SEASON, x20=X20):
    return (T - 1.0) * x20


X10 = on_manifold_x10()
_W, Y_JUNCTION, Z_JUNCTION = junction_constants()
ANALYTIC_T2 = T_SEASON - Y_JUNCTION
ANALYTIC_TERMINAL_RATIO = Y_JUNCTION + Z_JUNCTION  # x2(T)/x1(T)


def plant_model(t, x, u):
    veg, repro = x[0], x[1]
    alloc = u[0]
    dx = cs.vertcat(alloc * veg, (1.0 - alloc) * veg)
    y = Namespace(veg=veg, repro=repro, alloc=alloc)
    return dx, y


def build_ocp(T=T_SEASON, x10=None, x20=X20):
    if x10 is None:
        x10 = on_manifold_x10(T=T, x20=x20)
    sys = eco_opt.System(states=2, controls=1, model=plant_model, name="king_rg")
    ocp = eco_opt.Ocp()
    ocp.max([cs.log(sys.y.repro)])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "<=", sys.u[0], "<=", 1.0),
            (0.0, "<=", sys.y.veg),
            (1e-8, "<=", sys.y.repro),
            (0.0, "==", eco_opt.t_0(sys.t)),
            (T, "==", eco_opt.t_f(sys.t)),
            (x10, "==", eco_opt.t_0(sys.y.veg)),
            (x20, "==", eco_opt.t_0(sys.y.repro)),
        ],
    )
    return ocp, sys


def solve(
    T=T_SEASON,
    x10=None,
    x20=X20,
    control_intervals=120,
    polynomial_degree=4,
):
    if x10 is None:
        x10 = on_manifold_x10(T=T, x20=x20)
    ocp, sys = build_ocp(T=T, x10=x10, x20=x20)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: x10, sys.x[1]: x20, sys.u[0]: 0.6},
        ipopt_options={"max_iter": 3000, "tol": 1e-8},
        tf_guess=T,
    )
    return sol, sys


def analytic_t2(T=T_SEASON):
    _, y, _ = junction_constants()
    return T - y


def analytic_control(t, T=T_SEASON):
    """Closed-form u*(t) for a start-on-manifold trajectory."""
    t = np.asarray(t, dtype=float)
    t2 = analytic_t2(T)
    u = np.zeros_like(t, dtype=float)
    singular = t < t2
    u[singular] = 1.0 - 1.0 / (T - 1.0 - t[singular])
    return u


def analytic_trajectories(t, T=T_SEASON, x20=X20):
    """Closed-form states and control; assumes start-on-manifold ICs."""
    t = np.asarray(t, dtype=float)
    t2 = analytic_t2(T)
    x2_sing = x20 * np.exp(np.minimum(t, t2))
    veg = np.empty_like(t, dtype=float)
    repro = np.empty_like(t, dtype=float)
    singular = t <= t2
    repro[singular] = x20 * np.exp(t[singular])
    veg[singular] = (T - 1.0 - t[singular]) * repro[singular]
    x2_t2 = x20 * np.exp(t2)
    x1_t2 = (T - 1.0 - t2) * x2_t2
    repro[~singular] = x1_t2 * (t[~singular] - t2) + x2_t2
    veg[~singular] = x1_t2
    return veg, repro, analytic_control(t, T=T)


def analytic_objective(T=T_SEASON, x20=X20):
    """J* = int_0^T ln(x2) dt on the start-on-manifold solution."""
    w, y, _ = junction_constants()
    t2 = T - y
    # On [0, t2]: x2 = x20 e^t. On [t2, T]: x2 = x20 e^{t2} (1 + w (t-t2)).
    # int ln(1 + w τ) dτ from 0 to y equals w^2 (see README).
    return T * math.log(x20) + 0.5 * t2 * t2 + y * t2 + w * w


ANALYTIC_OBJECTIVE = analytic_objective()


def singular_ratio(t, veg, repro, T=T_SEASON):
    """x1 / ((T-1-t) x2); equals 1 on the singular arc.

    Returns NaN where the denominator is non-positive (t >= T-1).
    """
    t = np.asarray(t, dtype=float)
    denom = (T - 1.0 - t) * np.asarray(repro, dtype=float)
    out = np.full_like(t, np.nan, dtype=float)
    ok = denom > 1e-12
    out[ok] = np.asarray(veg, dtype=float)[ok] / denom[ok]
    return out


def t2_from_control(sol, u_thresh=0.2):
    """First interval-start at which u has dropped to the reproductive bang.

    After t2 the analytic ratio x1/((T-1-t)x2) leaves 1, then returns through
    1 as the denominator T-1-t shrinks, so a last-time-near-1 test is the
    wrong detector. The control drop is the junction.
    """
    u = np.asarray(sol.u[:, 0])
    t = np.asarray(sol.time_intervals[:-1])
    low = u < u_thresh
    if not np.any(low):
        return float(t[-1])
    return float(t[int(np.nonzero(low)[0][0])])


def t2_from_state_ratio(sol, T=T_SEASON, rel_tol=0.05):
    """Last time the singular ratio is ~1 *before* it first leaves that band.

    Equivalent to t2_from_control on a clean bang-singular-bang solution;
    kept as a state-only diagnostic.
    """
    t = np.asarray(sol.time)
    veg = np.asarray(sol.x[:, 0])
    repro = np.asarray(sol.x[:, 1])
    mask = t < (T - 1.0 - 1e-6)
    ratio = singular_ratio(t[mask], veg[mask], repro[mask], T=T)
    on_arc = np.abs(ratio - 1.0) < rel_tol
    if not np.any(on_arc):
        return float(t[0])
    # first departure after the opening on-arc run
    leaves = np.where(~on_arc)[0]
    if leaves.size == 0:
        return float(t[mask][-1])
    first_leave = int(leaves[0])
    if first_leave == 0:
        return float(t[mask][0])
    return float(t[mask][first_leave - 1])


if __name__ == "__main__":
    w, y, z = junction_constants()
    print("junction w, y, z =", w, y, z)
    print("analytic t2 =", ANALYTIC_T2, "terminal ratio =", ANALYTIC_TERMINAL_RATIO)
    print("analytic J* =", ANALYTIC_OBJECTIVE)
    sol, sys = solve()
    print(sol)
    print("OCP t2 (control) =", t2_from_control(sol),
          "(state-ratio) =", t2_from_state_ratio(sol))
    print("OCP x2(T)/x1(T) =", sol.x[-1, 1] / sol.x[-1, 0])
    print("OCP J =", sol.objective)
