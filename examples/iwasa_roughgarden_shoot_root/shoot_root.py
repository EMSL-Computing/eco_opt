"""Iwasa & Roughgarden (1984) shoot/root/fruit allocation, Liebig case.

Qualitative theorem (Iwasa & Roughgarden, TPB 25:78–105)
--------------------------------------------------------
An annual plant with two vegetative organs (shoot S, root R) and a
reproductive sink (fruit F) that maximises season-end fruit has a
three-phase optimum:

1. Unbalanced bang: put every photosynthate unit into the currently
   limiting organ until the shoot/root ratio reaches the balanced path.
2. Balanced-growth singular arc: S and R grow together at the ratio that
   maximises daily photosynthesis for the current vegetative biomass;
   nothing is sent to fruit.
3. Reproductive bang: a simultaneous stop of shoot and root growth, with
   all subsequent photosynthate going to fruit.

The original paper allows a general increasing concave production
function G(S, R). This example implements the **Liebig / minimum**
special case

    G(S, R) = min(α S, β R)

for which the balanced path, the balanced allocations, and the fruit
switch time are closed form. The three-phase *structure* is
Iwasa–Roughgarden; the production function is the tractable instance,
not a verbatim transcription of their most general G.

Problem::

    maximize   F(T)
    subject to dS/dt = u_S * G
               dR/dt = u_R * G
               dF/dt = (1 - u_S - u_R) * G
               G = min(α S, β R)
               0 <= u_S,  0 <= u_R,  u_S + u_R <= 1
               S(0)=S0, R(0)=R0, F(0)=0,  t in [0, T]

Closed form (Liebig)
--------------------
Balanced ratio:     R/S = α/β
Balanced controls:  u_S = β/(α+β),  u_R = α/(α+β)
Balanced RGR:       γ = α β / (α+β)
Fruit switch:       t_s = T - 1/γ     (Cohen rule with RGR = γ)

If the seedling is root-limited (β R0 < α S0), phase 1 is all-to-root
(u_R=1, S frozen, R(t) = R0 exp(β t)) until β R = α S0:

    t_bal = (1/β) ln(α S0 / (β R0))

Then S(t) = S0 exp(γ (t - t_bal)) on the balanced arc. After t_s the
vegetative organs freeze and fruit accumulates at the constant rate
G(t_s), so

    F(T) = G(t_s) / γ = α S0 exp(γ (t_s - t_bal)) / γ

Default demo: α=β=1 so γ=1/2, T=6, S0=1, R0=0.2 (root-deficient)
⇒ t_bal = ln(5) ≈ 1.609, t_s = 4, F(T) = 2 e^2 / √5.

The Liebig kink sits *on* the balanced singular arc, so a raw
``cs.fmin`` stalls IPOPT (exact Hessian undefined on a positive-length
interval). The dynamics use a standard hyperbolic smooth-min

    G_eps = 0.5 * (a + b - sqrt((a-b)^2 + eps^2)),   a=αS, b=βR

with ``eps = 1e-3``. As ``eps → 0`` this is ``min(a,b)``. Tests compare
against the *true* Liebig closed form, not against a smoothed analytic.
"""

import math

import casadi as cs
import numpy as np

import eco_opt
from eco_opt import Namespace

ALPHA = 1.0
BETA = 1.0
T_SEASON = 6.0
S0 = 1.0
R0 = 0.2
F0 = 0.0
LIEBIG_EPS = 1e-3  # hyperbolic smooth-min of min(αS, βR); see liebig_G


def balanced_rgr(alpha=ALPHA, beta=BETA):
    return alpha * beta / (alpha + beta)


def balanced_controls(alpha=ALPHA, beta=BETA):
    denom = alpha + beta
    return beta / denom, alpha / denom  # u_S, u_R


def fruit_switch_time(T=T_SEASON, alpha=ALPHA, beta=BETA):
    return T - 1.0 / balanced_rgr(alpha, beta)


def balance_time(S0=S0, R0=R0, alpha=ALPHA, beta=BETA):
    """Time at which a Liebig seedling first hits R/S = α/β.

    Returns 0 if the initial ratio is already balanced (within 1e-12).
    """
    shoot_lim = alpha * S0
    root_lim = beta * R0
    if abs(shoot_lim - root_lim) <= 1e-12:
        return 0.0
    if root_lim < shoot_lim:
        # root-deficient: all to root, R grows at rate β, S frozen
        return (1.0 / beta) * math.log((alpha * S0) / (beta * R0))
    # shoot-deficient: all to shoot, S grows at rate α, R frozen
    return (1.0 / alpha) * math.log((beta * R0) / (alpha * S0))


def analytic_terminal_fruit(
    T=T_SEASON, S0=S0, R0=R0, alpha=ALPHA, beta=BETA
):
    gamma = balanced_rgr(alpha, beta)
    t_s = fruit_switch_time(T, alpha, beta)
    t_bal = balance_time(S0, R0, alpha, beta)
    if t_bal >= t_s:
        raise ValueError(
            f"balance time {t_bal} is not before fruit switch {t_s}; "
            "the three-phase Liebig solution does not apply"
        )
    shoot_lim = alpha * S0
    root_lim = beta * R0
    if root_lim < shoot_lim - 1e-12:
        s_bal = S0
    elif shoot_lim < root_lim - 1e-12:
        s_bal = (beta / alpha) * R0
    else:
        s_bal = S0
    g_ts = alpha * s_bal * math.exp(gamma * (t_s - t_bal))
    return g_ts * (T - t_s)


GAMMA = balanced_rgr()
ANALYTIC_T_BAL = balance_time()
ANALYTIC_T_S = fruit_switch_time()
ANALYTIC_FRUIT = analytic_terminal_fruit()
U_S_BAL, U_R_BAL = balanced_controls()


def liebig_G(S, R, alpha=ALPHA, beta=BETA, eps=LIEBIG_EPS):
    """Hyperbolic smooth-min of min(αS, βR).

    ``0.5*(a+b-sqrt((a-b)^2+eps^2)) → min(a,b)`` as ``eps → 0``. The
    kink of the true Liebig function lies on the balanced singular arc,
    so the C^1 surrogate is the transcription that IPOPT can sit on.
    """
    a = alpha * S
    b = beta * R
    return 0.5 * (a + b - cs.sqrt((a - b) ** 2 + eps ** 2))


def shoot_root_model(t, x, u, alpha=ALPHA, beta=BETA, eps=LIEBIG_EPS):
    S, R, F = x[0], x[1], x[2]
    u_S, u_R = u[0], u[1]
    G = liebig_G(S, R, alpha=alpha, beta=beta, eps=eps)
    u_F = 1.0 - u_S - u_R
    y = Namespace(S=S, R=R, F=F, u_S=u_S, u_R=u_R, u_F=u_F, G=G)
    return [u_S * G, u_R * G, u_F * G], y


def build_ocp(
    T=T_SEASON,
    S0=S0,
    R0=R0,
    F0=F0,
    alpha=ALPHA,
    beta=BETA,
    eps=LIEBIG_EPS,
):
    def model(t, x, u):
        return shoot_root_model(t, x, u, alpha=alpha, beta=beta, eps=eps)

    sys = eco_opt.System(states=3, controls=2, model=model, name="iwasa_rg")
    ocp = eco_opt.Ocp()
    ocp.max(eco_opt.t_f(sys.y.F))
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "<=", sys.u[0], "<=", 1.0),
            (0.0, "<=", sys.u[1], "<=", 1.0),
            (0.0, "<=", sys.u[0] + sys.u[1], "<=", 1.0),
            (0.0, "<=", sys.y.S),
            (0.0, "<=", sys.y.R),
            (0.0, "<=", sys.y.F),
            (0.0, "==", eco_opt.t_0(sys.t)),
            (T, "==", eco_opt.t_f(sys.t)),
            (S0, "==", eco_opt.t_0(sys.y.S)),
            (R0, "==", eco_opt.t_0(sys.y.R)),
            (F0, "==", eco_opt.t_0(sys.y.F)),
        ],
    )
    return ocp, sys


def solve(
    T=T_SEASON,
    S0=S0,
    R0=R0,
    F0=F0,
    alpha=ALPHA,
    beta=BETA,
    eps=LIEBIG_EPS,
    control_intervals=80,
    polynomial_degree=3,
):
    ocp, sys = build_ocp(
        T=T, S0=S0, R0=R0, F0=F0, alpha=alpha, beta=beta, eps=eps
    )
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="radau",
        initial_guess={
            sys.x[0]: S0,
            sys.x[1]: R0,
            sys.x[2]: 0.0,
            sys.u[0]: 0.35,
            sys.u[1]: 0.45,
        },
        ipopt_options={
            "max_iter": 5000,
            "tol": 1e-6,
            "acceptable_tol": 1e-4,
            "acceptable_iter": 15,
            "hessian_approximation": "limited-memory",
        },
        tf_guess=T,
    )
    return sol, sys


def analytic_trajectories(
    t, T=T_SEASON, S0=S0, R0=R0, alpha=ALPHA, beta=BETA
):
    """Piecewise closed-form S, R, F, u_S, u_R on an arbitrary grid."""
    t = np.asarray(t, dtype=float)
    gamma = balanced_rgr(alpha, beta)
    t_s = fruit_switch_time(T, alpha, beta)
    t_bal = balance_time(S0, R0, alpha, beta)
    u_S_b, u_R_b = balanced_controls(alpha, beta)

    shoot_lim = alpha * S0
    root_lim = beta * R0
    root_deficient = root_lim < shoot_lim - 1e-12
    shoot_deficient = shoot_lim < root_lim - 1e-12

    S = np.empty_like(t, dtype=float)
    R = np.empty_like(t, dtype=float)
    F = np.empty_like(t, dtype=float)
    u_S = np.empty_like(t, dtype=float)
    u_R = np.empty_like(t, dtype=float)

    if root_deficient:
        s_bal, r_bal = S0, (alpha / beta) * S0
    elif shoot_deficient:
        s_bal, r_bal = (beta / alpha) * R0, R0
    else:
        s_bal, r_bal = S0, R0

    s_ts = s_bal * math.exp(gamma * (t_s - t_bal))
    r_ts = r_bal * math.exp(gamma * (t_s - t_bal))
    g_ts = alpha * s_ts

    for i, ti in enumerate(np.atleast_1d(t).astype(float)):
        if ti < t_bal:
            if root_deficient:
                S[i] = S0
                R[i] = R0 * math.exp(beta * ti)
                u_S[i], u_R[i] = 0.0, 1.0
            elif shoot_deficient:
                S[i] = S0 * math.exp(alpha * ti)
                R[i] = R0
                u_S[i], u_R[i] = 1.0, 0.0
            else:
                S[i], R[i] = S0, R0
                u_S[i], u_R[i] = u_S_b, u_R_b
            F[i] = 0.0
        elif ti < t_s:
            dt = ti - t_bal
            S[i] = s_bal * math.exp(gamma * dt)
            R[i] = r_bal * math.exp(gamma * dt)
            F[i] = 0.0
            u_S[i], u_R[i] = u_S_b, u_R_b
        else:
            S[i], R[i] = s_ts, r_ts
            F[i] = g_ts * (ti - t_s)
            u_S[i], u_R[i] = 0.0, 0.0
    return S, R, F, u_S, u_R


def fruit_switch_from_ocp(sol, thresh=0.4):
    """First interval-start at which vegetative allocation drops below thresh."""
    veg = np.asarray(sol.u[:, 0]) + np.asarray(sol.u[:, 1])
    t = np.asarray(sol.time_intervals[:-1])
    low = veg < thresh
    if not np.any(low):
        return float(t[-1])
    return float(t[int(np.nonzero(low)[0][0])])


def balance_time_from_ocp(sol, alpha=ALPHA, beta=BETA, rel_tol=0.06):
    """First coarse-grid time at which R/S is within rel_tol of α/β."""
    S = np.asarray(sol.x[:, 0])
    R = np.asarray(sol.x[:, 1])
    t = np.asarray(sol.time)
    target = alpha / beta
    ratio = R / np.maximum(S, 1e-12)
    on = np.abs(ratio / target - 1.0) < rel_tol
    if not np.any(on):
        return float(t[-1])
    return float(t[int(np.nonzero(on)[0][0])])


if __name__ == "__main__":
    print("gamma =", GAMMA)
    print("analytic t_bal, t_s, F(T) =", ANALYTIC_T_BAL, ANALYTIC_T_S, ANALYTIC_FRUIT)
    sol, sys = solve()
    print(sol)
    print("OCP t_bal =", balance_time_from_ocp(sol), "OCP t_s =", fruit_switch_from_ocp(sol))
    print("OCP F(T) =", sol.x[-1, 2], "S(T), R(T) =", sol.x[-1, 0], sol.x[-1, 1])
    print("OCP R(T)/S(T) =", sol.x[-1, 1] / sol.x[-1, 0], "target", ALPHA / BETA)
