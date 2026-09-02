"""Optimal-control solution results.

Provides the objective value, final time, solver status, named-signal
trajectories, and ``.signal(expr)`` sampling on the solution grid.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np
import casadi as cs

__all__ = ["OcpResults"]


class OcpResults:
    """Holds a solved NLP plus enough symbolic bookkeeping to reconstruct
    state/control trajectories and evaluate arbitrary signal expressions
    (``sys.y.rocket.velocity`` etc.) on the resulting grid.
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
        all_syms = []
        if npar:
            all_syms.append(Pk_sym)
        all_syms += [T0_sym, Tf_sym]
        all_syms.append(X_all[0])
        for k in range(N):
            all_syms.append(U_all[k])
            for j in range(d):
                all_syms.append(Xc_all[k][j])
            all_syms.append(X_all[k + 1])

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
        for k in range(N):
            u_vals.append(vals[idx]); idx += 1
            xj = []
            for j in range(d):
                xj.append(vals[idx]); idx += 1
            xc_vals.append(xj)
            x_vals.append(vals[idx]); idx += 1

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
        fn = cs.Function("sig", [self.sys.t, self.sys.x, self.sys.u], [expr])
        out = []
        for k in range(self.N):
            u_k = self.u[k] if self.u.size else np.zeros(self.sys.nu)
            out.append(float(fn(self.time[k], self.x[k], u_k)))
        # append terminal value using the last control (piecewise-constant hold)
        u_last = self.u[-1] if self.u.size else np.zeros(self.sys.nu)
        out.append(float(fn(self.time[-1], self.x[-1], u_last)))
        return np.array(out)

    def __repr__(self):
        return (
            f"OcpResults(success={self.success}, status={self.return_status!r}, "
            f"objective={self.objective!r}, t0={self.t0!r}, tf={self.tf!r})"
        )
