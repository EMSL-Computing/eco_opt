"""Optimal-control solution results.

Provides the objective value, final time, solver status, named-signal
trajectories, and ``.signal(expr)`` sampling on the solution grid.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np
import casadi as cs

__all__ = ["OcpResults", "PhaseTrajectory"]


class PhaseTrajectory:
    """State/control trajectories on one phase of a (possibly multi-phase) OCP."""

    def __init__(self, sys, t0, tf, time, x, u, time_intervals, N, d, tau=None, z=None):
        self.sys = sys
        self.t0 = float(t0)
        self.tf = float(tf)
        self.time = np.asarray(time)
        self.x = np.asarray(x)
        self.u = np.asarray(u)
        nz = getattr(sys, "nz", 0)
        n_t = self.x.shape[0] if self.x.ndim else 0
        self.z = np.asarray(z) if z is not None else np.zeros((n_t, nz))
        self.time_intervals = np.asarray(time_intervals)
        self.N = int(N)
        self.d = int(d)
        self.tau = np.asarray(tau) if tau is not None else np.zeros(0)

    def signal(self, expr):
        fn = cs.Function("sig", [self.sys.t, self.sys.x, self.sys.u, self.sys.z], [expr])
        out = []
        for k in range(self.N):
            u_k = self.u[k] if self.u.size else np.zeros(self.sys.nu)
            z_k = self.z[k] if self.z.size else np.zeros(self.sys.nz)
            out.append(float(fn(self.time[k], self.x[k], u_k, z_k)))
        u_last = self.u[-1] if self.u.size else np.zeros(self.sys.nu)
        z_last = self.z[-1] if self.z.size else np.zeros(self.sys.nz)
        out.append(float(fn(self.time[-1], self.x[-1], u_last, z_last)))
        return np.array(out)

    def __repr__(self):
        return (
            f"PhaseTrajectory(sys={self.sys.name!r}, t0={self.t0!r}, "
            f"tf={self.tf!r}, N={self.N})"
        )


class OcpResults:
    """Holds a solved NLP plus enough symbolic bookkeeping to reconstruct
    state/control trajectories and evaluate arbitrary signal expressions
    (``sys.y.rocket.velocity`` etc.) on the resulting grid.

    Single-phase solutions also carry shadow prices recovered from the NLP
    multipliers of the collocation equations:

    * ``costate`` (N*d, nx): costate ``lam = dJ/dx`` at the collocation
      points, in the sense of the declared objective. For ``max``, a
      positive value is the gain in ``J`` per extra unit of that state.
    * ``costate_time`` (N*d,): the collocation-point times.
    * ``hamiltonian`` (N*d,): ``H = L + lam' f`` at the same points; for an
      autonomous problem it is constant along an optimum, and zero when
      the final time is free.

    Multi-phase results set these to ``None``.
    """

    def __init__(
        self,
        sys,
        sol,
        stats,
        N: int,
        d: int,
        npar: int,
        tau,
        T0_sym,
        Tf_sym,
        Pk_sym,
        X_all: List,
        U_all: List,
        Xc_all: List[List],
        w,
        objective_sense: Optional[str],
        Z_all=None,
        Zc_all=None,
        costate=None,
        hamiltonian_fn=None,
    ):
        self.sys = sys
        self.stats = stats
        self.N = N
        self.d = d
        self.tau = np.asarray(tau)

        w_opt = np.array(sol["x"]).ravel()

        # Rebuild a numeric-lookup function from the symbolic decomposition:
        # extract each symbol's slice from the flat w vector via a CasADi
        # Function (robust against however variables were interleaved).
        nz = int(getattr(sys, "nz", 0))
        use_z = bool(nz and Z_all)
        all_syms = []
        if npar:
            all_syms.append(Pk_sym)
        all_syms += [T0_sym, Tf_sym]
        all_syms.append(X_all[0])
        for k in range(N):
            all_syms.append(U_all[k])
            if use_z:
                all_syms.append(Z_all[k])
            for j in range(d):
                all_syms.append(Xc_all[k][j])
                if use_z and Zc_all:
                    all_syms.append(Zc_all[k][j])
            all_syms.append(X_all[k + 1])
        if use_z:
            all_syms.append(Z_all[N])

        extract_fn = cs.Function("extract", [w], all_syms)
        vals = extract_fn(w_opt)
        if not isinstance(vals, (list, tuple)):
            vals = [vals]
        vals = [np.array(v).ravel() for v in vals]

        idx = 0
        if npar:
            p_val = vals[idx]
            idx += 1
        else:
            p_val = np.zeros(0)
        t0_val = float(vals[idx][0]); idx += 1
        tf_val = float(vals[idx][0]); idx += 1
        x0_val = vals[idx]; idx += 1

        u_vals = []
        xc_vals = []  # xc_vals[k][j]
        x_vals = [x0_val]
        z_vals = []
        zc_vals = []
        for k in range(N):
            u_vals.append(vals[idx]); idx += 1
            if use_z:
                z_vals.append(vals[idx]); idx += 1
            xj = []
            zj = []
            for j in range(d):
                xj.append(vals[idx]); idx += 1
                if use_z and Zc_all:
                    zj.append(vals[idx]); idx += 1
            xc_vals.append(xj)
            zc_vals.append(zj)
            x_vals.append(vals[idx]); idx += 1
        if use_z:
            z_vals.append(vals[idx]); idx += 1

        self.p = p_val
        self.t0 = t0_val
        self.tf = tf_val
        self.objective_sense = objective_sense
        f_val = float(sol["f"])
        self.objective = -f_val if objective_sense == "max" else f_val

        h = (tf_val - t0_val) / N

        # Build full time grid: interval starts + collocation points + interval ends
        t_nodes = []
        x_nodes = []
        for k in range(N):
            tk = t0_val + k * h
            t_nodes.append(tk)
            x_nodes.append(x_vals[k])
            for j in range(d):
                t_nodes.append(tk + self.tau[j + 1] * h)
                x_nodes.append(xc_vals[k][j])
        t_nodes.append(tf_val)
        x_nodes.append(x_vals[N])

        self.time_nodes = np.array(t_nodes)
        self.x_nodes = np.array(x_nodes)  # shape (num_nodes, nx)

        # Piecewise-constant control grid (u constant on each interval)
        self.time_intervals = np.array([t0_val + k * h for k in range(N)] + [tf_val])
        self.u_intervals = np.array(u_vals) if u_vals else np.zeros((0, sys.nu))

        # State at the coarse (interval-boundary) grid -- useful for direct
        # comparison against reference solutions sampled once per interval.
        self.time = np.array([t0_val + k * h for k in range(N)] + [tf_val])
        self.x = np.array(x_vals)  # shape (N+1, nx)
        self.u = self.u_intervals  # shape (N, nu), piecewise-constant
        if use_z:
            self.z = np.array(z_vals)  # shape (N+1, nz) at interval boundaries
        else:
            self.z = np.zeros((N + 1, 0))

        # Shadow prices (costates) and Hamiltonian at the collocation points.
        self.costate_time = None
        self.costate = None
        self.hamiltonian = None
        if costate is not None:
            t_c = [t0_val + (k + self.tau[j + 1]) * h for k in range(N) for j in range(d)]
            self.costate_time = np.array(t_c)
            self.costate = np.asarray(costate)
            if hamiltonian_fn is not None:
                H = []
                n = 0
                for k in range(N):
                    for j in range(d):
                        z_kj = zc_vals[k][j] if (use_z and zc_vals[k]) else np.zeros(nz)
                        H.append(float(hamiltonian_fn(
                            t_c[n], xc_vals[k][j], u_vals[k], z_kj, p_val, self.costate[n]
                        )))
                        n += 1
                self.hamiltonian = np.array(H)
        self.phases = [self]
        self.switch_times = np.zeros(0)

    @classmethod
    def from_phases(cls, phase_trajs, stats, objective, objective_sense):
        """Assemble a multi-phase result. ``.phases`` holds each stage;
        ``.t0``/``.tf``/``.objective`` are the overall problem; ``.time``
        and ``.x`` are concatenated when every phase has the same ``nx``.
        """
        obj = cls.__new__(cls)
        obj.stats = stats
        obj.objective_sense = objective_sense
        obj.objective = float(objective)
        obj.phases = list(phase_trajs)
        obj.t0 = float(phase_trajs[0].t0)
        obj.tf = float(phase_trajs[-1].tf)
        obj.sys = phase_trajs[-1].sys
        obj.N = int(sum(p.N for p in phase_trajs))
        obj.d = int(phase_trajs[0].d)
        obj.tau = phase_trajs[0].tau
        obj.p = np.zeros(0)
        obj.switch_times = np.array([p.tf for p in phase_trajs[:-1]], dtype=float)
        same_nx = len({p.x.shape[1] if p.x.ndim == 2 else 0 for p in phase_trajs}) == 1
        if same_nx:
            times = [phase_trajs[0].time]
            xs = [phase_trajs[0].x]
            us = [phase_trajs[0].u]
            tints = [phase_trajs[0].time_intervals]
            for p in phase_trajs[1:]:
                times.append(p.time[1:])
                xs.append(p.x[1:])
                us.append(p.u)
                tints.append(p.time_intervals[1:])
            obj.time = np.concatenate(times)
            obj.x = np.concatenate(xs, axis=0)
            obj.u = np.concatenate(us, axis=0) if us[0].ndim == 2 else np.zeros((0, 0))
            obj.time_intervals = np.concatenate(tints)
            zs = [phase_trajs[0].z]
            for p in phase_trajs[1:]:
                zs.append(p.z[1:] if getattr(p, "z", None) is not None and p.z.size else p.z)
            try:
                obj.z = np.concatenate(zs, axis=0)
            except ValueError:
                obj.z = phase_trajs[-1].z
        else:
            last = phase_trajs[-1]
            obj.time = last.time
            obj.x = last.x
            obj.u = last.u
            obj.time_intervals = last.time_intervals
            obj.z = getattr(last, "z", np.zeros((last.x.shape[0], 0)))
        obj.time_nodes = obj.time
        obj.x_nodes = obj.x
        obj.costate_time = None
        obj.costate = None
        obj.hamiltonian = None
        obj.u_intervals = obj.u
        return obj

    # -- status -------------------------------------------------------------
    @property
    def success(self) -> bool:
        return bool(self.stats.get("success", False))

    @property
    def return_status(self) -> str:
        return self.stats.get("return_status", "unknown")

    # -- signal sampling ------------------------------------------------------
    def state(self, index: int) -> np.ndarray:
        """Return the coarse-grid trajectory of state component `index`."""
        return self.x[:, index]

    def control(self, index: int = 0) -> np.ndarray:
        """Return the piecewise-constant control trajectory for component `index`."""
        return self.u[:, index] if self.u.size else np.zeros(0)

    def signal(self, expr_or_index, kind: str = "state") -> np.ndarray:
        """Retrieve a named/derived signal sampled on the coarse time grid.

        Pass a symbolic expression built from ``sys.x``/``sys.u`` (e.g.
        ``sys.y.rocket.velocity``) and get a numeric array aligned with
        ``self.time`` (states) or ``self.time_intervals[:-1]`` (controls).
        """
        if isinstance(expr_or_index, int):
            return self.state(expr_or_index) if kind == "state" else self.control(expr_or_index)

        expr = expr_or_index
        fn = cs.Function("sig", [self.sys.t, self.sys.x, self.sys.u, self.sys.z], [expr])
        out = []
        z_traj = getattr(self, "z", np.zeros((self.N + 1, getattr(self.sys, "nz", 0))))
        for k in range(self.N):
            u_k = self.u[k] if self.u.size else np.zeros(self.sys.nu)
            z_k = z_traj[k] if z_traj.size else np.zeros(self.sys.nz)
            out.append(float(fn(self.time[k], self.x[k], u_k, z_k)))
        u_last = self.u[-1] if self.u.size else np.zeros(self.sys.nu)
        z_last = z_traj[-1] if z_traj.size else np.zeros(self.sys.nz)
        out.append(float(fn(self.time[-1], self.x[-1], u_last, z_last)))
        return np.array(out)

    def __repr__(self):
        return (
            f"OcpResults(success={self.success}, status={self.return_status!r}, "
            f"objective={self.objective!r}, t0={self.t0!r}, tf={self.tf!r})"
        )
