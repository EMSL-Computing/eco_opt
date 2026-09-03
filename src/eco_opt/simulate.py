"""Open-loop replay of a control and shrinking/receding-horizon closed loop.

``simulate`` integrates ``ẋ = f(t,x,u)`` (CasADi CVODES / IDAS) with a
given input. That is the plant. ``closed_loop`` repeatedly solves the OCP
from the current plant state, applies the first piece of ``u*`` for one
sampling interval, and integrates the plant — concatenation of open-loop
prefixes. If the plant is the OCP model and nothing is disturbed,
shrinking-horizon closed loop reprints the open-loop tape.
"""

from __future__ import annotations

from typing import Callable, Optional

import numpy as np
import casadi as cs

from .constraints import classify_constraint
from .expressions import t_0, t_f
from .results import OcpResults

__all__ = [
    "SimulationResults",
    "ClosedLoopResults",
    "simulate",
    "closed_loop",
    "control_from_solution",
]


class SimulationResults:
    """Plant trajectory: states on ``time``, controls constant on each interval."""

    def __init__(self, sys, time, x, u, z=None, p=None):
        self.sys = sys
        self.time = np.asarray(time, dtype=float)
        self.x = np.asarray(x, dtype=float)
        self.u = np.asarray(u, dtype=float)
        n_t = self.time.size
        nz = getattr(sys, "nz", 0)
        self.z = (
            np.asarray(z, dtype=float)
            if z is not None
            else np.zeros((n_t, nz))
        )
        self.p = None if p is None else np.asarray(p, dtype=float).ravel()
        self.t0 = float(self.time[0])
        self.tf = float(self.time[-1])

    def __repr__(self):
        return (
            f"SimulationResults(sys={self.sys.name!r}, t0={self.t0!r}, "
            f"tf={self.tf!r}, n={self.time.size})"
        )


class ClosedLoopResults:
    """Plant trajectory under receding/shrinking-horizon application of ``u*``."""

    def __init__(
        self,
        sys,
        time,
        x,
        u,
        z=None,
        p=None,
        statuses=None,
        objectives=None,
        horizon="shrinking",
    ):
        self.sys = sys
        self.time = np.asarray(time, dtype=float)
        self.x = np.asarray(x, dtype=float)
        self.u = np.asarray(u, dtype=float)
        n_t = self.time.size
        nz = getattr(sys, "nz", 0)
        self.z = (
            np.asarray(z, dtype=float)
            if z is not None
            else np.zeros((n_t, nz))
        )
        self.p = None if p is None else np.asarray(p, dtype=float).ravel()
        self.statuses = list(statuses or [])
        self.objectives = np.asarray(objectives, dtype=float) if objectives is not None else np.zeros(0)
        self.horizon = horizon
        self.t0 = float(self.time[0])
        self.tf = float(self.time[-1])
        self.success = all(
            s in ("Solve_Succeeded", "Solved_To_Acceptable_Level")
            for s in self.statuses
        ) if self.statuses else False

    def __repr__(self):
        return (
            f"ClosedLoopResults(horizon={self.horizon!r}, success={self.success}, "
            f"steps={len(self.statuses)}, t0={self.t0!r}, tf={self.tf!r})"
        )


def control_from_solution(sol: OcpResults) -> Callable[[float], np.ndarray]:
    """Piecewise-constant ``u(t)`` from a collocation solution."""
    edges = np.asarray(sol.time_intervals, dtype=float)
    uu = np.asarray(sol.u, dtype=float)
    nu = sol.sys.nu
    if uu.size == 0:
        z = np.zeros(nu)

        def u_of_t(_t):
            return z.copy()

        return u_of_t

    def u_of_t(t):
        t = float(t)
        if t >= edges[-1]:
            k = uu.shape[0] - 1
        elif t <= edges[0]:
            k = 0
        else:
            k = int(np.searchsorted(edges, t, side="right") - 1)
            k = min(max(k, 0), uu.shape[0] - 1)
        return np.asarray(uu[k], dtype=float).ravel()

    return u_of_t


def _as_u_fn(u, nu: int) -> Callable[[float], np.ndarray]:
    if u is None:
        z = np.zeros(nu)

        def u_of_t(_t):
            return z.copy()

        return u_of_t
    if isinstance(u, OcpResults):
        return control_from_solution(u)
    if isinstance(u, ClosedLoopResults):
        edges = u.time
        uu = u.u

        def u_of_t(t):
            t = float(t)
            if t >= edges[-1]:
                k = uu.shape[0] - 1
            elif t <= edges[0]:
                k = 0
            else:
                k = int(np.searchsorted(edges, t, side="right") - 1)
                k = min(max(k, 0), uu.shape[0] - 1)
            return np.asarray(uu[k], dtype=float).ravel()

        return u_of_t
    if callable(u):
        def u_of_t(t):
            val = np.asarray(u(float(t)), dtype=float).ravel()
            if val.size == 1 and nu > 1:
                val = np.repeat(val, nu)
            if val.size != nu:
                raise ValueError(f"u(t) returned {val.size} entries; expected {nu}")
            return val

        return u_of_t
    arr = np.asarray(u, dtype=float).ravel()
    if arr.size == 1 and nu > 1:
        arr = np.repeat(arr, nu)
    if arr.size != nu:
        raise ValueError(f"constant u has {arr.size} entries; expected {nu}")

    def u_of_t(_t, _arr=arr):
        return _arr.copy()

    return u_of_t


def _pack_p(sys, u_val, p_val):
    parts = []
    if sys.nu:
        parts.append(np.asarray(u_val, dtype=float).ravel()[: sys.nu])
    if sys.np:
        if p_val is None:
            raise ValueError("system has parameters; pass p=...")
        parts.append(np.asarray(p_val, dtype=float).ravel()[: sys.np])
    if not parts:
        return np.zeros(0)
    return np.concatenate(parts)


def _dae_symbols(sys):
    p_parts = []
    if sys.nu:
        p_parts.append(sys.u)
    if sys.np:
        p_parts.append(sys.p)
    p_sym = cs.vertcat(*p_parts) if p_parts else cs.MX.sym("p_empty", 0)
    dae = {"x": sys.x, "ode": sys.dx, "t": sys.t, "p": p_sym}
    if sys.nz:
        dae["z"] = sys.z
        dae["alg"] = sys.dz
    return dae


def _consistent_z(sys, t, x, u_val, p_val, z0):
    if not sys.nz:
        return np.zeros(0)
    if z0 is not None:
        return np.asarray(z0, dtype=float).ravel()[: sys.nz]
    z = cs.MX.sym("z0", sys.nz)
    g = cs.substitute(sys.dz, sys.t, float(t))
    g = cs.substitute(g, sys.x, np.asarray(x, dtype=float).ravel())
    if sys.nu:
        g = cs.substitute(g, sys.u, np.asarray(u_val, dtype=float).ravel())
    if sys.np:
        g = cs.substitute(g, sys.p, np.asarray(p_val, dtype=float).ravel())
    g = cs.substitute(g, sys.z, z)
    rf = cs.rootfinder("z0", "newton", {"x": z, "g": g})
    return np.array(rf(np.zeros(sys.nz))).ravel()


def _step_integrator(sys, t0, tf):
    dae = _dae_symbols(sys)
    plugin = "idas" if sys.nz else "cvodes"
    opts = {"abstol": 1e-8, "reltol": 1e-6}
    return cs.integrator("plant", plugin, dae, float(t0), float(tf), opts)


def simulate(
    sys,
    x0,
    t_grid,
    u=None,
    p=None,
    z0=None,
) -> SimulationResults:
    """Integrate the system ODE/DAE on ``t_grid`` with input ``u``.

    ``u`` may be a constant, a callable ``u(t)``, or an :class:`OcpResults`
    (piecewise-constant collocation control). Time in ``f`` is the
    integrator clock, so interpolants of exogenous inputs keep working.
    """
    t_grid = np.asarray(t_grid, dtype=float).ravel()
    if t_grid.size < 2:
        raise ValueError("t_grid must have at least two points")
    x0 = np.asarray(x0, dtype=float).ravel()
    if x0.size != sys.nx:
        raise ValueError(f"x0 has {x0.size} entries; system has {sys.nx} states")
    u_fn = _as_u_fn(u, sys.nu)
    xs = [x0.copy()]
    us = []
    zs = []
    xk = x0.copy()
    uk0 = u_fn(t_grid[0])
    zk = _consistent_z(sys, t_grid[0], xk, uk0, p, z0)
    zs.append(zk.copy())
    for k in range(t_grid.size - 1):
        t0, tf = float(t_grid[k]), float(t_grid[k + 1])
        uk = u_fn(t0)
        us.append(uk)
        p_num = _pack_p(sys, uk, p)
        step = _step_integrator(sys, t0, tf)
        arg = {"x0": xk, "p": p_num}
        if sys.nz:
            arg["z0"] = zk
        out = step(**arg)
        xk = np.array(out["xf"]).ravel()
        xs.append(xk.copy())
        if sys.nz:
            zk = np.array(out["zf"]).ravel()
        zs.append(zk.copy())
    u_arr = np.vstack(us) if us else np.zeros((0, sys.nu))
    return SimulationResults(sys, t_grid, np.vstack(xs), u_arr, z=np.vstack(zs), p=p)


def _clock_from_constraints(sys, constraints):
    t0 = 0.0
    tf = None
    for c in constraints:
        pc = classify_constraint(sys, c)
        if pc.slot and pc.slot[0] == "t":
            if pc.kind == "box_initial" and np.isfinite(pc.lb) and abs(pc.lb - pc.ub) < 1e-14:
                t0 = float(pc.lb)
            elif pc.kind == "box_final" and np.isfinite(pc.lb) and abs(pc.lb - pc.ub) < 1e-14:
                tf = float(pc.lb)
            elif pc.kind == "box":
                if np.isfinite(pc.lb):
                    t0 = float(pc.lb)
                if np.isfinite(pc.ub):
                    tf = float(pc.ub)
    return t0, tf


def _retune_clock_and_state(sys, constraints, t0, tf, x0):
    """Copy ``constraints``, replacing t0/tf pins and initial-state equalities."""
    x0 = np.asarray(x0, dtype=float).ravel()
    new = []
    saw_t0 = saw_tf = False
    saw_x = [False] * sys.nx
    for c in constraints:
        pc = classify_constraint(sys, c)
        if pc.kind == "box_initial" and pc.slot and pc.slot[0] == "t":
            new.append((float(t0), "==", t_0(sys.t)))
            saw_t0 = True
            continue
        if pc.kind == "box_final" and pc.slot and pc.slot[0] == "t":
            new.append((float(tf), "==", t_f(sys.t)))
            saw_tf = True
            continue
        if pc.kind == "box_initial" and pc.slot and pc.slot[0] == "x":
            i = pc.slot[1]
            new.append((float(x0[i]), "==", t_0(sys.x[i])))
            saw_x[i] = True
            continue
        new.append(c)
    if not saw_t0:
        new.append((float(t0), "==", t_0(sys.t)))
    if not saw_tf:
        new.append((float(tf), "==", t_f(sys.t)))
    for i in range(sys.nx):
        if not saw_x[i]:
            new.append((float(x0[i]), "==", t_0(sys.x[i])))
    return new


def _state_path_box(sys, constraints):
    lo = np.full(sys.nx, -np.inf)
    hi = np.full(sys.nx, np.inf)
    for c in constraints:
        pc = classify_constraint(sys, c)
        if pc.kind == "box" and pc.slot and pc.slot[0] == "x":
            i = pc.slot[1]
            lo[i] = max(lo[i], pc.lb)
            hi[i] = min(hi[i], pc.ub)
    return lo, hi


def _clip_state_to_box(x, lo, hi):
    x = np.asarray(x, dtype=float).ravel().copy()
    for i in range(x.size):
        if np.isfinite(lo[i]) and np.isfinite(hi[i]) and lo[i] > hi[i]:
            continue
        if np.isfinite(lo[i]):
            x[i] = max(x[i], lo[i])
        if np.isfinite(hi[i]):
            x[i] = min(x[i], hi[i])
    return x


def _paint_previous(compiled, t0, tf, u_fn, x_now, x_end):
    N, d = compiled.N, compiled.d
    h = (float(tf) - float(t0)) / N
    tau = np.asarray(compiled.tau, dtype=float)
    nu = compiled.sys.nu
    nx = compiled.sys.nx
    if nu and u_fn is not None:
        for i in range(nu):
            for k, idx in enumerate(compiled.idx.u[i]):
                tmid = t0 + (k + 0.5) * h
                compiled.w0[idx] = float(np.asarray(u_fn(tmid)).ravel()[i])
    x_now = np.asarray(x_now, dtype=float).ravel()
    x_end = np.asarray(x_end, dtype=float).ravel() if x_end is not None else x_now
    span = max(tf - t0, 1e-16)
    for i in range(nx):
        slots = compiled.idx.x[i]
        ptr = 0
        compiled.w0[slots[ptr]] = float(x_now[i])
        ptr += 1
        for k in range(N):
            t_left = t0 + k * h
            for j in range(d):
                t = t_left + h * tau[j + 1]
                a = (t - t0) / span
                compiled.w0[slots[ptr]] = float((1 - a) * x_now[i] + a * x_end[i])
                ptr += 1
            t_right = t0 + (k + 1) * h
            a = (t_right - t0) / span
            compiled.w0[slots[ptr]] = float((1 - a) * x_now[i] + a * x_end[i])
            ptr += 1


def closed_loop(
    ocp,
    x0,
    dt: float,
    t_final: Optional[float] = None,
    *,
    horizon: str = "shrinking",
    prediction_horizon: Optional[float] = None,
    plant=None,
    p=None,
    z0=None,
    control_intervals: Optional[int] = None,
    polynomial_degree: int = 3,
    collocation_points: str = "radau",
    ipopt_options: Optional[dict] = None,
    initial_guess: Optional[dict] = None,
) -> ClosedLoopResults:
    """Apply the OCP in a shrinking or receding horizon on a plant integrator.

    At each sampling instant ``t_k`` the compiled (or rebuilt) NLP is solved
    from the current plant state. Only the first piecewise-constant control
    is applied for ``dt``; the plant is integrated with that hold.

    ``horizon='shrinking'`` keeps the original terminal time (time-consistent
    with one finite-horizon OCP). ``horizon='receding'`` always looks
    ``prediction_horizon`` ahead (a different controller).

    Multi-phase problems are not supported.
    """
    if ocp._phases:
        raise NotImplementedError("closed_loop is single-phase")
    if ocp.system is None:
        raise RuntimeError("Ocp.st(...) must be called before closed_loop()")
    if horizon not in ("shrinking", "receding"):
        raise ValueError("horizon must be 'shrinking' or 'receding'")
    sys = ocp.system
    plant = plant if plant is not None else sys
    if plant.nx != sys.nx or plant.nu != sys.nu:
        raise ValueError("plant must have the same (nx, nu) as the OCP system")
    dt = float(dt)
    if dt <= 0:
        raise ValueError("dt must be positive")
    x0 = np.asarray(x0, dtype=float).ravel()
    t_start, tf_con = _clock_from_constraints(sys, ocp.constraints)
    if t_final is None:
        t_final = tf_con
    if t_final is None:
        raise ValueError("t_final is not pinned in the OCP; pass t_final=...")
    t_final = float(t_final)
    if horizon == "receding":
        if prediction_horizon is None:
            prediction_horizon = t_final - t_start
        prediction_horizon = float(prediction_horizon)
        if prediction_horizon <= 0:
            raise ValueError("prediction_horizon must be positive")

    ipopt_options = dict(ipopt_options or {})
    ipopt_options.setdefault("print_level", 0)
    ipopt_options.setdefault("max_iter", 2000)

    t_k = float(t_start)
    x_lo, x_hi = _state_path_box(sys, ocp.constraints)
    x_k = _clip_state_to_box(x0, x_lo, x_hi)
    z_k = z0
    times = [t_k]
    xs = [x_k.copy()]
    us = []
    zs_hist = []
    statuses = []
    objectives = []
    u_fn_prev = None
    x_end_prev = None
    base_constraints = list(ocp.constraints)

    n_steps = int(np.round((t_final - t_start) / dt))
    if n_steps < 1:
        raise ValueError("t_final - t0 must be at least one dt")

    try:
        for _step in range(n_steps):
            remaining = t_final - t_k
            if remaining <= 1e-12:
                break
            step_dt = min(dt, remaining)
            if horizon == "shrinking":
                tf_ocp = t_final
                n_int = int(max(1, round((tf_ocp - t_k) / dt)))
            else:
                tf_ocp = t_k + prediction_horizon
                n_int = (
                    int(control_intervals)
                    if control_intervals is not None
                    else int(max(1, round(prediction_horizon / dt)))
                )
            if control_intervals is not None and horizon == "shrinking":
                n_int = int(control_intervals)

            ocp.st(
                system=sys,
                constraints=_retune_clock_and_state(
                    sys, base_constraints, t_k, tf_ocp, x_k
                ),
            )
            guess = dict(initial_guess or {})
            for i in range(sys.nx):
                guess.setdefault(sys.x[i], float(x_k[i]))
            ocp.build(
                control_intervals=n_int,
                polynomial_degree=polynomial_degree,
                collocation_points=collocation_points,
                initial_guess=guess if guess else None,
                tf_guess=tf_ocp,
            )
            if u_fn_prev is not None:
                _paint_previous(ocp._compiled, t_k, tf_ocp, u_fn_prev, x_k, x_end_prev)
            sol = ocp.optimize(
                initial_guess=None,
                ipopt_options=ipopt_options,
                tf_guess=tf_ocp,
                warm_start=False,
            )
            statuses.append(sol.return_status)
            objectives.append(float(sol.objective))
            if not sol.success and sol.return_status != "Solved_To_Acceptable_Level":
                break
            u_apply = (
                np.asarray(sol.u[0], dtype=float).ravel()
                if sol.u.size
                else np.zeros(sys.nu)
            )
            us.append(u_apply)
            u_fn_prev = control_from_solution(sol)
            x_end_prev = np.asarray(sol.x[-1], dtype=float).ravel()
            p_use = p
            if p_use is None and getattr(sol, "p", None) is not None and np.size(sol.p):
                p_use = sol.p
            z_use = z_k
            if z_use is None and getattr(sol, "z", None) is not None and np.size(sol.z):
                z_use = sol.z[0]
            sim = simulate(
                plant,
                x_k,
                [t_k, t_k + step_dt],
                u=u_apply,
                p=p_use,
                z0=z_use,
            )
            if plant.nz and not zs_hist:
                zs_hist.append(np.asarray(sim.z[0], dtype=float).ravel())
            x_k = _clip_state_to_box(sim.x[-1], x_lo, x_hi)
            z_k = sim.z[-1] if sim.z.size else z_k
            t_k = t_k + step_dt
            times.append(t_k)
            xs.append(x_k.copy())
            if plant.nz:
                zs_hist.append(np.asarray(sim.z[-1], dtype=float).ravel())
    finally:
        ocp.st(system=sys, constraints=base_constraints)

    z_arr = (
        np.vstack(zs_hist)
        if zs_hist
        else np.zeros((len(times), plant.nz))
    )
    u_arr = np.vstack(us) if us else np.zeros((0, sys.nu))
    return ClosedLoopResults(
        plant,
        times,
        np.vstack(xs),
        u_arr,
        z=z_arr,
        p=p,
        statuses=statuses,
        objectives=objectives,
        horizon=horizon,
    )
