"""Optimal control problem definition and direct-collocation transcription.

Declare an objective with ``min``/``max``, constraints with ``st()``
(including ``t_0``/``t_f`` timepoint wrapping and box/path/initial/final
classification), and solve by Radau or Legendre orthogonal collocation.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Union

import numpy as np
import casadi as cs

from .expressions import TimepointExpression
from .constraints import classify_constraint, ParsedConstraint
from .collocation import CollocationCoefficients
from .results import OcpResults, PhaseTrajectory
from .phase import Phase

__all__ = ["Ocp"]

_MIN_PHASE_DURATION = 1e-6


def _depends_on_system(expr, sys) -> bool:
    raw = expr.expr if isinstance(expr, TimepointExpression) else expr
    try:
        if sys.nx and cs.depends_on(raw, sys.x):
            return True
        if sys.nu and cs.depends_on(raw, sys.u):
            return True
        if sys.nz and cs.depends_on(raw, sys.z):
            return True
        if cs.depends_on(raw, sys.t):
            return True
    except Exception:
        return False
    return False


def _knots_from_tf_guess(tf_guess, n_phases: int) -> np.ndarray:
    if isinstance(tf_guess, (int, float)):
        return np.linspace(0.0, float(tf_guess), n_phases + 1)
    arr = np.asarray(tf_guess, dtype=float).ravel()
    if arr.size == n_phases + 1:
        return arr
    if arr.size == n_phases:
        return np.concatenate([[0.0], np.cumsum(arr)])
    raise ValueError(
        "tf_guess for a multi-phase problem must be a scalar overall tf, "
        f"a sequence of {n_phases + 1} knot times, or {n_phases} durations"
    )


def _expr_equal(a, b) -> bool:
    try:
        return bool(cs.is_equal(a, b, 2))
    except Exception:
        return False


def _mx_len(expr) -> int:
    return int(expr.numel())


def _clip_guess(w0, lbw, ubw):
    lo = np.where(np.isfinite(lbw), lbw, w0)
    hi = np.where(np.isfinite(ubw), ubw, w0)
    return np.minimum(np.maximum(w0, lo), hi)


def _solver_opts(ipopt_options):
    opts = {"ipopt": dict(ipopt_options or {})}
    opts["ipopt"].setdefault("print_level", 0)
    opts.setdefault("print_time", 0)
    return opts


class _BoundIndex:
    """Maps box slots and constraint lists onto flat NLP bound vectors."""

    def __init__(self, nx, nu, npar, nz=0):
        self.p = [[] for _ in range(npar)]
        self.t0 = None
        self.tf = None
        self.x0 = [None] * nx
        self.xf = [None] * nx
        self.x = [[] for _ in range(nx)]
        self.u = [[] for _ in range(nu)]
        self.z0 = [None] * nz
        self.zf = [None] * nz
        self.z = [[] for _ in range(nz)]
        self.g_initial = []
        self.g_path = []
        self.g_final = []


class _CompiledNlp:
    """Compiled single-phase collocation NLP: change numeric bounds, re-solve."""

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.solver = None
        self.solver_key = None
        self.w_last = None

    def write_bounds(self):
        lbw, ubw = self.lbw, self.ubw
        box, b0, bf = self.box, self.box_initial, self.box_final
        idx = self.idx
        sys = self.sys
        for i in range(sys.np):
            lo, hi = box[("p", i)]
            for k in idx.p[i]:
                lbw[k] = lo
                ubw[k] = hi
        lbw[idx.t0] = self.t0_bounds[0]
        ubw[idx.t0] = self.t0_bounds[1]
        lbw[idx.tf] = self.tf_bounds[0]
        ubw[idx.tf] = self.tf_bounds[1]
        for i in range(sys.nx):
            lo, hi = box[("x", i)]
            for k in idx.x[i]:
                lbw[k] = lo
                ubw[k] = hi
            k0 = idx.x0[i]
            lbw[k0] = max(lo, b0[("x", i)][0])
            ubw[k0] = min(hi, b0[("x", i)][1])
            kf = idx.xf[i]
            if kf is not None:
                lbw[kf] = max(lo, bf[("x", i)][0])
                ubw[kf] = min(hi, bf[("x", i)][1])
        for i in range(sys.nu):
            lo, hi = box[("u", i)]
            for k in idx.u[i]:
                lbw[k] = lo
                ubw[k] = hi
        for i in range(sys.nz):
            lo, hi = box[("z", i)]
            for k in idx.z[i]:
                lbw[k] = lo
                ubw[k] = hi
            if idx.z0[i] is not None:
                lbw[idx.z0[i]] = max(lo, b0.get(("z", i), (-np.inf, np.inf))[0])
                ubw[idx.z0[i]] = min(hi, b0.get(("z", i), (-np.inf, np.inf))[1])
            if idx.zf[i] is not None:
                lbw[idx.zf[i]] = max(lo, bf.get(("z", i), (-np.inf, np.inf))[0])
                ubw[idx.zf[i]] = min(hi, bf.get(("z", i), (-np.inf, np.inf))[1])
        for n, pc in enumerate(self.initial_constraints):
            for k in idx.g_initial[n]:
                self.lbg[k] = pc.lb
                self.ubg[k] = pc.ub
        for n, pc in enumerate(self.path_constraints):
            for k in idx.g_path[n]:
                self.lbg[k] = pc.lb
                self.ubg[k] = pc.ub
        for n, pc in enumerate(self.final_constraints):
            for k in idx.g_final[n]:
                self.lbg[k] = pc.lb
                self.ubg[k] = pc.ub

    def apply_constant_guess(self, x_guess, u_guess, p_guess, tf_guess, z_guess=None):
        w0 = self.w0
        idx = self.idx
        if z_guess is None:
            z_guess = getattr(self, "z_guess", [])
        for i, val in enumerate(p_guess):
            for k in idx.p[i]:
                w0[k] = val
        w0[idx.t0] = 0.0
        w0[idx.tf] = float(tf_guess)
        for i, val in enumerate(x_guess):
            w0[idx.x0[i]] = val
            for k in idx.x[i]:
                w0[k] = val
        for i, val in enumerate(u_guess):
            for k in idx.u[i]:
                w0[k] = val
        for i, val in enumerate(z_guess):
            for k in idx.z[i]:
                w0[k] = val
        self.x_guess = list(x_guess)
        self.u_guess = list(u_guess)
        self.p_guess = list(p_guess)
        self.z_guess = list(z_guess)
        self.tf_guess = float(tf_guess)

    def set_constraints(self, constraints):
        sys = self.sys
        parsed = [classify_constraint(sys, c) for c in constraints]
        for pc in parsed:
            if pc.kind == "box":
                if pc.slot[0] == "t":
                    self.t0_bounds = [pc.lb, pc.ub]
                    self.tf_bounds = [pc.lb, pc.ub]
                else:
                    self.box[pc.slot] = (pc.lb, pc.ub)
                    self._guess_from_equality(pc.slot, pc.lb, pc.ub, at="path")
            elif pc.kind == "box_initial":
                if pc.slot[0] == "t":
                    self.t0_bounds = [pc.lb, pc.ub]
                else:
                    self.box_initial[pc.slot] = (pc.lb, pc.ub)
                    self._guess_from_equality(pc.slot, pc.lb, pc.ub, at="initial")
            elif pc.kind == "box_final":
                if pc.slot[0] == "t":
                    self.tf_bounds = [pc.lb, pc.ub]
                else:
                    self.box_final[pc.slot] = (pc.lb, pc.ub)
            elif pc.kind in ("path", "initial", "final"):
                bucket = {
                    "path": self.path_constraints,
                    "initial": self.initial_constraints,
                    "final": self.final_constraints,
                }[pc.kind]
                matched = False
                for existing in bucket:
                    if _expr_equal(existing.expr, pc.expr):
                        existing.lb = pc.lb
                        existing.ub = pc.ub
                        matched = True
                if not matched:
                    raise ValueError(
                        "constraint is not in the compiled NLP; add it in "
                        "st() and call build() again"
                    )
            else:
                raise ValueError(f"unsupported constraint kind {pc.kind!r}")
        self.write_bounds()

    def _guess_from_equality(self, slot, lb, ub, at):
        if not (np.isfinite(lb) and np.isfinite(ub) and abs(ub - lb) < 1e-14):
            return
        val = float(lb)
        kind, i = slot
        if kind == "p":
            self.p_guess[i] = val
            for k in self.idx.p[i]:
                self.w0[k] = val
        elif kind == "x" and at == "initial":
            self.x_guess[i] = val
            self.w0[self.idx.x0[i]] = val
        elif kind == "u":
            self.u_guess[i] = val
            for k in self.idx.u[i]:
                self.w0[k] = val
        elif kind == "z":
            self.z_guess[i] = val
            for k in self.idx.z[i]:
                self.w0[k] = val

    def ensure_solver(self, ipopt_options):
        key = tuple(sorted((ipopt_options or {}).items()))
        if self.solver is None or key != self.solver_key:
            self.solver = cs.nlpsol("solver", "ipopt", self.nlp, _solver_opts(ipopt_options))
            self.solver_key = key
        return self.solver

    def run(self, ipopt_options, warm_start, initial_guess, tf_guess, sense):
        if initial_guess is not None:
            self._apply_guess_dict(initial_guess)
            warm_start = False
        if tf_guess is not None:
            self.w0[self.idx.tf] = float(tf_guess)
            self.tf_guess = float(tf_guess)
        if warm_start and self.w_last is not None:
            self.w0 = self.w_last.copy()
            if tf_guess is not None:
                self.w0[self.idx.tf] = float(tf_guess)
        self.w0 = _clip_guess(self.w0, self.lbw, self.ubw)
        solver = self.ensure_solver(ipopt_options)
        sol = solver(
            x0=self.w0, lbx=self.lbw, ubx=self.ubw, lbg=self.lbg, ubg=self.ubg
        )
        stats = solver.stats()
        self.w_last = np.array(sol["x"]).ravel()
        return OcpResults(
            sys=self.sys,
            sol=sol,
            stats=stats,
            N=self.N,
            d=self.d,
            npar=self.sys.np,
            tau=self.tau,
            T0_sym=self.T0,
            Tf_sym=self.Tf,
            Pk_sym=self.Pk,
            X_all=self.X_all,
            U_all=self.U_all,
            Xc_all=self.Xc_all,
            Z_all=getattr(self, "Z_all", None),
            Zc_all=getattr(self, "Zc_all", None),
            w=self.w_sym,
            objective_sense=sense,
        )

    def _apply_guess_dict(self, initial_guess):
        sys = self.sys
        for key, val in initial_guess.items():
            xi = sys.state_index(key)
            ui = sys.control_index(key)
            pi = sys.parameter_index(key)
            zi = sys.algebraic_index(key)
            if xi is not None:
                self.x_guess[xi] = float(val)
            elif ui is not None:
                self.u_guess[ui] = float(val)
            elif pi is not None:
                self.p_guess[pi] = float(val)
            elif zi is not None:
                self.z_guess[zi] = float(val)
            elif _expr_equal(sys.t, key):
                pass
            else:
                raise ValueError(
                    f"initial_guess key {key!r} is not a component of system {sys!r}"
                )
        self.apply_constant_guess(
            self.x_guess, self.u_guess, self.p_guess, self.tf_guess, self.z_guess
        )


class _CompiledMultiphase:
    """Compiled sequential-phase NLP: re-run IPOPT without re-discretizing."""

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.solver = None
        self.solver_key = None
        self.w_last = None

    def set_constraints(self, constraints):
        raise NotImplementedError(
            "set_constraints is single-phase; change Phase constraints and call build() again"
        )

    def ensure_solver(self, ipopt_options):
        key = tuple(sorted((ipopt_options or {}).items()))
        if self.solver is None or key != self.solver_key:
            self.solver = cs.nlpsol(
                "solver", "ipopt", self.nlp, _solver_opts(ipopt_options)
            )
            self.solver_key = key
        return self.solver

    def _apply_tf_guess(self, tf_guess):
        knots = _knots_from_tf_guess(tf_guess, self.n_phases)
        for i in range(self.n_phases):
            self.w0[self.t0_idx[i]] = float(knots[i])
            self.w0[self.tf_idx[i]] = float(knots[i + 1])

    def _apply_guess_dict(self, initial_guess):
        for ph, sys in enumerate(self.systems):
            for key, val in initial_guess.items():
                xi = sys.state_index(key)
                ui = sys.control_index(key)
                if xi is not None:
                    for k in self.x_idx[ph][xi]:
                        self.w0[k] = float(val)
                elif ui is not None:
                    for k in self.u_idx[ph][ui]:
                        self.w0[k] = float(val)

    def run(self, ipopt_options, warm_start, initial_guess, tf_guess, sense):
        if initial_guess is not None:
            self._apply_guess_dict(initial_guess)
            warm_start = False
        if tf_guess is not None:
            self._apply_tf_guess(tf_guess)
        if warm_start and self.w_last is not None:
            self.w0 = self.w_last.copy()
            if tf_guess is not None:
                self._apply_tf_guess(tf_guess)
        self.w0 = _clip_guess(self.w0, self.lbw, self.ubw)
        solver = self.ensure_solver(ipopt_options)
        sol = solver(
            x0=self.w0, lbx=self.lbw, ubx=self.ubw, lbg=self.lbg, ubg=self.ubg
        )
        stats = solver.stats()
        self.w_last = np.array(sol["x"]).ravel()
        return self._results(sol, stats, sense)

    def _results(self, sol, stats, sense):
        w_opt = np.array(sol["x"]).ravel()
        trajs = []
        cc = self.cc
        w_sym = self.w_sym
        for pack in self.phase_pack:
            all_syms = [pack["T0"], pack["Tf"], pack["X_all"][0]]
            for k in range(pack["N"]):
                all_syms.append(pack["U_all"][k])
                for j in range(pack["d"]):
                    all_syms.append(pack["Xc_all"][k][j])
                all_syms.append(pack["X_all"][k + 1])
            extract_fn = cs.Function("extract", [w_sym], all_syms)
            vals = extract_fn(w_opt)
            if not isinstance(vals, (list, tuple)):
                vals = [vals]
            vals = [np.array(v).ravel() for v in vals]
            t0_val = float(vals[0][0])
            tf_val = float(vals[1][0])
            idx = 2
            x_vals = [vals[idx]]
            idx += 1
            u_vals = []
            Nloc = pack["N"]
            dloc = pack["d"]
            for k in range(Nloc):
                u_vals.append(vals[idx])
                idx += 1
                idx += dloc
                x_vals.append(vals[idx])
                idx += 1
            t_coarse = np.array(
                [t0_val + k * (tf_val - t0_val) / Nloc for k in range(Nloc)]
                + [tf_val]
            )
            u_arr = np.array(u_vals) if pack["nu"] else np.zeros((Nloc, 0))
            trajs.append(
                PhaseTrajectory(
                    sys=pack["sys"],
                    t0=t0_val,
                    tf=tf_val,
                    time=t_coarse,
                    x=np.array(x_vals),
                    u=u_arr,
                    time_intervals=t_coarse,
                    N=Nloc,
                    d=dloc,
                    tau=cc.tau,
                )
            )
        f_val = float(sol["f"])
        objective = -f_val if sense == "max" else f_val
        return OcpResults.from_phases(trajs, stats, objective, sense)


class Ocp:
    """Optimal control problem: objective + constraints + solve.

    Usage::

        ocp = Ocp()
        ocp.max(t_f(sig.height))
        ocp.st(
            system=sys,
            constraints=[
                (0, '==', t_0(time)),
                (0, '==', t_0(sig.velocity)),
                (0, '<=', sig.velocity, '<=', float('inf')),
            ],
        )
        sol = ocp.solve(control_intervals=59, polynomial_degree=5,
                        collocation_points='legendre')
        sim = ocp.simulate(sol)          # integrate ẋ=f(t,x,u*)
        cl = ocp.closed_loop(x0, dt, t_final=T)  # shrinking horizon

    Bound / parameter sweeps without rebuilding the collocation graph::

        ocp.build(control_intervals=40)
        sol_a = ocp.optimize()
        ocp.set_constraints([(new_p, "==", sys.p[0])])
        sol_b = ocp.optimize()  # warm-starts from sol_a

    Multi-phase (sequential stages, optional continuity of time and state)::

        p1 = Phase(sys_a, constraints=[...])
        p2 = Phase(sys_b, constraints=[...])
        ocp.st(phases=[p1, p2], continuity=True)
        sol = ocp.solve(tf_guess=5.0)
        # sol.phases[i], sol.switch_times
    """

    def __init__(self):
        self._sense: Optional[str] = None  # 'min' or 'max'
        self._objective_terms: List = []
        self.system = None
        self.constraints: List = []
        self._phases: Optional[List[Phase]] = None
        self._continuity = True
        self._compiled: Optional[_CompiledNlp] = None

    def _invalidate(self):
        self._compiled = None

    # -- problem definition ----------------------------------------------
    def min(self, expression):
        self._invalidate()
        self._sense = "min"
        self._objective_terms = list(expression) if isinstance(expression, (list, tuple)) and not isinstance(expression, TimepointExpression) else [expression]

    def max(self, expression):
        self._invalidate()
        self._sense = "max"
        self._objective_terms = list(expression) if isinstance(expression, (list, tuple)) and not isinstance(expression, TimepointExpression) else [expression]

    def st(
        self,
        system=None,
        systems=None,
        constraints=None,
        phases=None,
        continuity: bool = True,
    ):
        self._invalidate()
        if phases is not None:
            phases = list(phases)
            if not phases:
                raise ValueError("st(phases=...) requires a non-empty sequence")
            for p in phases:
                if not isinstance(p, Phase):
                    raise TypeError(
                        "st(phases=...) expects eco_opt.Phase instances; "
                        "got {0!r}".format(type(p))
                    )
            self._phases = phases
            self._continuity = bool(continuity)
            self.system = phases[0].system
            self.constraints = []
            return
        self._phases = None
        self._continuity = True
        sys = system if system is not None else systems
        if isinstance(sys, (list, tuple)):
            raise ValueError(
                "multiple systems must be passed as phases=[Phase(sys, constraints), ...]"
            )
        self.system = sys
        if constraints is None:
            raise ValueError("st() requires constraints=[...]")
        self.constraints = list(constraints)

    # -- transcription + solve --------------------------------------------
    def solve(
        self,
        control_intervals: int = 100,
        polynomial_degree: int = 5,
        collocation_points: str = "legendre",
        initial_guess: Optional[dict] = None,
        ipopt_options: Optional[dict] = None,
        tf_guess: float = 1.0,
    ) -> OcpResults:
        """Build the collocation NLP and solve it (one-shot).

        For a bound or parameter sweep on a fixed mesh, call ``build`` once
        then ``set_constraints`` / ``optimize``.
        """
        self.build(
            control_intervals=control_intervals,
            polynomial_degree=polynomial_degree,
            collocation_points=collocation_points,
            initial_guess=initial_guess,
            tf_guess=tf_guess,
        )
        return self.optimize(
            initial_guess=initial_guess,
            ipopt_options=ipopt_options,
            tf_guess=tf_guess,
            warm_start=False,
        )

    def build(
        self,
        control_intervals: int = 100,
        polynomial_degree: int = 5,
        collocation_points: str = "legendre",
        initial_guess: Optional[dict] = None,
        tf_guess: float = 1.0,
    ):
        """Compile the collocation NLP. Mesh is frozen until the next ``build``.

        Numeric bounds and pinned parameters can then be changed with
        ``set_constraints`` without re-discretizing (single-phase).
        """
        if self._phases:
            self._compiled = self._compile_multiphase(
                control_intervals=control_intervals,
                polynomial_degree=polynomial_degree,
                collocation_points=collocation_points,
                initial_guess=initial_guess,
                tf_guess=tf_guess,
            )
            return self
        if self.system is None:
            raise RuntimeError("Ocp.st(...) must be called before build()")
        self._compiled = self._compile_singlephase(
            control_intervals=control_intervals,
            polynomial_degree=polynomial_degree,
            collocation_points=collocation_points,
            initial_guess=initial_guess,
            tf_guess=tf_guess,
        )
        return self

    def set_constraints(self, constraints):
        """Replace numeric bounds on constraints already declared in ``st()``.

        Same tuple syntax as ``st(constraints=...)``. Box bounds on states,
        controls, parameters, and time, and path/boundary expressions that
        were in the original problem, can be retuned. Adding a new path
        expression requires ``st()`` + ``build()``.
        """
        if self._compiled is None:
            raise RuntimeError("build() must be called before set_constraints()")
        self._compiled.set_constraints(list(constraints))
        return self

    def optimize(
        self,
        initial_guess: Optional[dict] = None,
        ipopt_options: Optional[dict] = None,
        tf_guess: Optional[float] = None,
        warm_start: bool = True,
    ) -> OcpResults:
        """Run IPOPT on the compiled NLP.

        After the first successful call, ``warm_start=True`` (default) uses
        the previous primal as ``x0``. Pass ``initial_guess`` to override
        with per-component constants.
        """
        if self._compiled is None:
            raise RuntimeError("build() must be called before optimize()")
        return self._compiled.run(
            ipopt_options=ipopt_options,
            warm_start=warm_start,
            initial_guess=initial_guess,
            tf_guess=tf_guess,
            sense=self._sense,
        )

    def simulate(self, sol, t_grid=None, x0=None, plant=None, p=None, z0=None):
        """Integrate the plant with a given control (open-loop replay).

        Default plant is this OCP's ``system``, default ``u`` is ``sol.u``
        (piecewise constant), default ``x0`` is ``sol.x[0]``, default grid
        is ``sol.time``.
        """
        from .simulate import simulate as _simulate

        sys = plant if plant is not None else self.system
        if t_grid is None:
            t_grid = sol.time
        if x0 is None:
            x0 = sol.x[0]
        if p is None and getattr(sol, "p", None) is not None and np.size(sol.p):
            p = sol.p
        if z0 is None and getattr(sol, "z", None) is not None and np.size(sol.z):
            z0 = sol.z[0]
        return _simulate(sys, x0, t_grid, u=sol, p=p, z0=z0)

    def closed_loop(self, x0, dt, t_final=None, **kwargs):
        """Shrinking- or receding-horizon application of this OCP on a plant.

        See ``eco_opt.closed_loop`` for arguments. Mutates and restores
        this object's constraints around the loop.
        """
        from .simulate import closed_loop as _closed_loop

        return _closed_loop(self, x0, dt, t_final=t_final, **kwargs)

    def _compile_singlephase(
        self,
        control_intervals: int,
        polynomial_degree: int,
        collocation_points: str,
        initial_guess: Optional[dict],
        tf_guess: float,
    ) -> _CompiledNlp:
        sys = self.system
        N = int(control_intervals)
        d = int(polynomial_degree)
        cc = CollocationCoefficients(d, collocation_points)
        C, D, B = cc.C, cc.D, cc.B

        nx, nu, nz, npar = sys.nx, sys.nu, sys.nz, sys.np
        if nz and sys.dz is None:
            raise ValueError("algebraic system is missing residual dz")

        f = cs.Function("f", [sys.t, sys.x, sys.u, sys.z, sys.p], [sys.dx])
        alg = (
            cs.Function("alg", [sys.t, sys.x, sys.u, sys.z, sys.p], [sys.dz])
            if nz
            else None
        )

        x_guess = [1.0] * nx
        u_guess = [1.0] * nu
        p_guess = [1.0] * npar
        z_guess = [1.0] * nz
        if initial_guess:
            for key, val in initial_guess.items():
                xi = sys.state_index(key)
                ui = sys.control_index(key)
                pi = sys.parameter_index(key)
                zi = sys.algebraic_index(key)
                if xi is not None:
                    x_guess[xi] = float(val)
                elif ui is not None:
                    u_guess[ui] = float(val)
                elif pi is not None:
                    p_guess[pi] = float(val)
                elif zi is not None:
                    z_guess[zi] = float(val)
                elif cs.is_equal(sys.t, key, 1):
                    pass
                else:
                    raise ValueError(f"initial_guess key {key!r} is not a component of system {sys!r}")

        parsed = [classify_constraint(sys, c) for c in self.constraints]

        box_initial = {("x", i): (-np.inf, np.inf) for i in range(nx)}
        box_final = {("x", i): (-np.inf, np.inf) for i in range(nx)}
        box = {("x", i): (-np.inf, np.inf) for i in range(nx)}
        box.update({("u", i): (-np.inf, np.inf) for i in range(nu)})
        box.update({("z", i): (-np.inf, np.inf) for i in range(nz)})
        box.update({("p", i): (-np.inf, np.inf) for i in range(npar)})
        box_initial.update({("z", i): (-np.inf, np.inf) for i in range(nz)})
        box_final.update({("z", i): (-np.inf, np.inf) for i in range(nz)})
        t0_bounds = [-np.inf, np.inf]
        tf_bounds = [-np.inf, np.inf]

        path_constraints: List[ParsedConstraint] = []
        initial_constraints: List[ParsedConstraint] = []
        final_constraints: List[ParsedConstraint] = []

        for pc in parsed:
            if pc.kind == "box":
                if pc.slot[0] == "t":
                    t0_bounds = [max(t0_bounds[0], pc.lb), min(t0_bounds[1], pc.ub)]
                    tf_bounds = [max(tf_bounds[0], pc.lb), min(tf_bounds[1], pc.ub)]
                else:
                    lo, hi = box[pc.slot]
                    box[pc.slot] = (max(lo, pc.lb), min(hi, pc.ub))
            elif pc.kind == "box_initial":
                if pc.slot[0] == "t":
                    t0_bounds = [max(t0_bounds[0], pc.lb), min(t0_bounds[1], pc.ub)]
                else:
                    lo, hi = box_initial.get(pc.slot, (-cs.inf, cs.inf))
                    box_initial[pc.slot] = (max(lo, pc.lb), min(hi, pc.ub))
            elif pc.kind == "box_final":
                if pc.slot[0] == "t":
                    tf_bounds = [max(tf_bounds[0], pc.lb), min(tf_bounds[1], pc.ub)]
                else:
                    lo, hi = box_final.get(pc.slot, (-cs.inf, cs.inf))
                    box_final[pc.slot] = (max(lo, pc.lb), min(hi, pc.ub))
            elif pc.kind == "path":
                path_constraints.append(pc)
            elif pc.kind == "initial":
                initial_constraints.append(pc)
            elif pc.kind == "final":
                final_constraints.append(pc)

        mayer_terms = []
        lagrange_terms = []
        for term in self._objective_terms:
            if isinstance(term, TimepointExpression):
                mayer_terms.append(term)
            else:
                lagrange_terms.append(term)
        has_lagrange = len(lagrange_terms) > 0
        if has_lagrange:
            lagrange_fn = cs.Function(
                "L", [sys.t, sys.x, sys.u, sys.z, sys.p], [sum(lagrange_terms)]
            )

        w = []
        w0 = []
        lbw = []
        ubw = []
        g = []
        lbg = []
        ubg = []
        idx = _BoundIndex(nx, nu, npar, nz)
        cursor = 0
        g_cursor = 0

        def push_w(lo, hi, guess):
            nonlocal cursor
            k = cursor
            w0.append(guess)
            lbw.append(lo)
            ubw.append(hi)
            cursor += 1
            return k

        def push_g(expr, lo, hi):
            nonlocal g_cursor
            n = _mx_len(expr)
            g.append(expr)
            if n == 1:
                lbg.append(lo)
                ubg.append(hi)
                sl = [g_cursor]
                g_cursor += 1
            else:
                lbg.extend([lo] * n)
                ubg.extend([hi] * n)
                sl = list(range(g_cursor, g_cursor + n))
                g_cursor += n
            return sl

        if npar:
            Pk = cs.MX.sym("P", npar)
            w.append(Pk)
            for i in range(npar):
                lo, hi = box[("p", i)]
                idx.p[i].append(push_w(lo, hi, p_guess[i]))
        else:
            Pk = cs.MX.sym("P", 0)

        T0 = cs.MX.sym("T0")
        Tf = cs.MX.sym("Tf")
        w += [T0, Tf]
        idx.t0 = push_w(t0_bounds[0], t0_bounds[1], 0.0)
        idx.tf = push_w(tf_bounds[0], tf_bounds[1], tf_guess)

        h_interval = (Tf - T0) / N

        Xk = cs.MX.sym("X0", nx)
        w.append(Xk)
        for i in range(nx):
            lo, hi = box[("x", i)]
            lo = max(lo, box_initial[("x", i)][0])
            hi = min(hi, box_initial[("x", i)][1])
            k = push_w(lo, hi, x_guess[i])
            idx.x0[i] = k
            idx.x[i].append(k)

        def eval_expr_at(expr, t_val, x_val, u_val, z_val):
            fn = cs.Function("tmp", [sys.t, sys.x, sys.u, sys.z, sys.p], [expr])
            return fn(t_val, x_val, u_val, z_val, Pk)

        def push_z(name, lo_hi_guess=None):
            """Add nz algebraic NLP variables; return the MX vector."""
            Z = cs.MX.sym(name, nz) if nz else cs.MX.sym(name, 0)
            if nz:
                w.append(Z)
                for i in range(nz):
                    lo, hi = box[("z", i)]
                    ki = push_w(lo, hi, z_guess[i])
                    idx.z[i].append(ki)
            return Z

        z_empty = cs.MX.zeros(0) if not nz else None

        idx.g_initial = []
        X_all = [Xk]
        U_all = []
        Xc_all = []
        Z_all = []
        Zc_all = []
        idx.g_path = [[] for _ in path_constraints]

        J = 0.0

        for k in range(N):
            Uk = cs.MX.sym(f"U_{k}", nu) if nu else cs.MX.sym(f"U_{k}", 0)
            w.append(Uk)
            for i in range(nu):
                lo, hi = box[("u", i)]
                idx.u[i].append(push_w(lo, hi, u_guess[i]))
            U_all.append(Uk)

            Zk = push_z(f"Z_{k}")
            if nz:
                if k == 0:
                    for i in range(nz):
                        idx.z0[i] = idx.z[i][-1]
                        lo, hi = box[("z", i)]
                        lo = max(lo, box_initial[("z", i)][0])
                        hi = min(hi, box_initial[("z", i)][1])
                        lbw[idx.z0[i]] = lo
                        ubw[idx.z0[i]] = hi
                tk = T0 + k * h_interval
                push_g(alg(tk, Xk, Uk, Zk, Pk), 0.0, 0.0)
            else:
                Zk = z_empty
            Z_all.append(Zk)

            if k == 0:
                for ic in initial_constraints:
                    val = eval_expr_at(ic.expr, T0, Xk, Uk, Zk)
                    idx.g_initial.append(push_g(val, ic.lb, ic.ub))

            Xkj = []
            Zkj = []
            for j in range(d):
                Xkj_ = cs.MX.sym(f"X_{k}_{j}", nx)
                w.append(Xkj_)
                for i in range(nx):
                    lo, hi = box[("x", i)]
                    idx.x[i].append(push_w(lo, hi, x_guess[i]))
                Xkj.append(Xkj_)
                Zkj_ = push_z(f"Zc_{k}_{j}")
                if not nz:
                    Zkj_ = z_empty
                Zkj.append(Zkj_)

            Xk_end = D[0] * Xk
            tk = T0 + k * h_interval
            # Collocation skips tau=0. Piecewise-constant Uk is defined at
            # the interval left end, so mixed g(x,u) must be applied there.
            for pi, pc in enumerate(path_constraints):
                val = eval_expr_at(pc.expr, tk, Xk, Uk, Zk)
                idx.g_path[pi].extend(push_g(val, pc.lb, pc.ub))
            for j in range(d):
                xp = C[0, j + 1] * Xk
                for r in range(d):
                    xp = xp + C[r + 1, j + 1] * Xkj[r]

                t_cj = tk + cc.tau[j + 1] * h_interval
                Zcj = Zkj[j] if nz else z_empty
                fj = f(t_cj, Xkj[j], Uk, Zcj, Pk)
                defect = h_interval * fj - xp
                push_g(defect, 0.0, 0.0)
                if nz:
                    push_g(alg(t_cj, Xkj[j], Uk, Zcj, Pk), 0.0, 0.0)

                Xk_end = Xk_end + D[j + 1] * Xkj[j]

                if has_lagrange:
                    J = J + B[j + 1] * h_interval * lagrange_fn(
                        t_cj, Xkj[j], Uk, Zcj, Pk
                    )

                for pi, pc in enumerate(path_constraints):
                    val = eval_expr_at(pc.expr, t_cj, Xkj[j], Uk, Zcj)
                    idx.g_path[pi].extend(push_g(val, pc.lb, pc.ub))

            Xk_next = cs.MX.sym(f"X_{k+1}", nx)
            w.append(Xk_next)
            is_last = k == N - 1
            for i in range(nx):
                lo, hi = box[("x", i)]
                if is_last:
                    lo = max(lo, box_final[("x", i)][0])
                    hi = min(hi, box_final[("x", i)][1])
                ki = push_w(lo, hi, x_guess[i])
                idx.x[i].append(ki)
                if is_last:
                    idx.xf[i] = ki

            push_g(Xk_end - Xk_next, 0.0, 0.0)

            Xk = Xk_next
            X_all.append(Xk)
            Xc_all.append(Xkj)
            Zc_all.append(Zkj)

        Zf = z_empty
        if nz:
            Zf = push_z("Zf")
            for i in range(nz):
                idx.zf[i] = idx.z[i][-1]
                lo, hi = box[("z", i)]
                lo = max(lo, box_final[("z", i)][0])
                hi = min(hi, box_final[("z", i)][1])
                lbw[idx.zf[i]] = lo
                ubw[idx.zf[i]] = hi
            u_last = U_all[-1]
            push_g(alg(Tf, Xk, u_last, Zf, Pk), 0.0, 0.0)
        Z_all.append(Zf)

        idx.g_final = []
        u_last = U_all[-1] if U_all else (cs.MX.zeros(nu) if nu else cs.MX.zeros(0))
        z_last = Zf if nz else z_empty
        for pi, pc in enumerate(path_constraints):
            val = eval_expr_at(pc.expr, Tf, Xk, u_last, z_last)
            idx.g_path[pi].extend(push_g(val, pc.lb, pc.ub))
        for fc in final_constraints:
            val = eval_expr_at(fc.expr, Tf, Xk, u_last, z_last)
            idx.g_final.append(push_g(val, fc.lb, fc.ub))

        for mt in mayer_terms:
            t_val = T0 if mt.point == "0" else Tf
            x_val = X_all[0] if mt.point == "0" else X_all[-1]
            u_dummy = U_all[0] if (mt.point == "0" and U_all) else (
                U_all[-1] if U_all else (cs.MX.zeros(nu) if nu else cs.MX.zeros(0))
            )
            z_dummy = Z_all[0] if (mt.point == "0" and nz) else z_last
            val = eval_expr_at(mt.expr, t_val, x_val, u_dummy, z_dummy)
            J = J + val

        if self._sense == "max":
            J = -J

        w_sym = cs.vertcat(*w)
        g_sym = cs.vertcat(*g) if g else cs.MX.sym("g_empty", 0)

        compiled = _CompiledNlp(
            sys=sys,
            N=N,
            d=d,
            tau=cc.tau,
            nlp={"f": J, "x": w_sym, "g": g_sym},
            w_sym=w_sym,
            T0=T0,
            Tf=Tf,
            Pk=Pk,
            X_all=X_all,
            U_all=U_all,
            Xc_all=Xc_all,
            Z_all=Z_all if nz else None,
            Zc_all=Zc_all if nz else None,
            lbw=np.asarray(lbw, dtype=float),
            ubw=np.asarray(ubw, dtype=float),
            lbg=np.asarray(lbg, dtype=float),
            ubg=np.asarray(ubg, dtype=float),
            w0=np.asarray(w0, dtype=float),
            box=box,
            box_initial=box_initial,
            box_final=box_final,
            t0_bounds=t0_bounds,
            tf_bounds=tf_bounds,
            path_constraints=path_constraints,
            initial_constraints=initial_constraints,
            final_constraints=final_constraints,
            idx=idx,
            x_guess=x_guess,
            u_guess=u_guess,
            p_guess=p_guess,
            z_guess=z_guess,
            tf_guess=float(tf_guess),
        )
        return compiled

    def _compile_multiphase(
        self,
        control_intervals: int,
        polynomial_degree: int,
        collocation_points: str,
        initial_guess: Optional[dict],
        tf_guess: float,
    ) -> _CompiledMultiphase:
        phases = self._phases
        n_ph = len(phases)
        d = int(polynomial_degree)
        cc = CollocationCoefficients(d, collocation_points)
        C, D, B = cc.C, cc.D, cc.B
        knots = _knots_from_tf_guess(tf_guess, n_ph)

        for p in phases:
            if p.system.nz:
                raise NotImplementedError(
                    "algebraic states are not supported in multi-phase transcription"
                )
            if p.system.np:
                raise NotImplementedError(
                    "parameters are not yet supported in multi-phase transcription"
                )

        # Split Mayer / Lagrange terms onto the phase whose system they use.
        mayer_by_phase = [[] for _ in range(n_ph)]
        lagrange_by_phase = [[] for _ in range(n_ph)]
        for term in self._objective_terms:
            assigned = False
            for i, p in enumerate(phases):
                if _depends_on_system(term, p.system):
                    if isinstance(term, TimepointExpression):
                        mayer_by_phase[i].append(term)
                    else:
                        lagrange_by_phase[i].append(term)
                    assigned = True
                    break
            if not assigned:
                if isinstance(term, TimepointExpression):
                    mayer_by_phase[-1].append(term)
                else:
                    lagrange_by_phase[-1].append(term)

        w, w0, lbw, ubw = [], [], [], []
        g, lbg, ubg = [], [], []
        J = 0.0
        Pk = cs.MX.sym("P", 0)
        cursor = 0
        t0_idx, tf_idx = [], []
        x_idx, u_idx = [], []
        systems = []

        def push_w(lo, hi, guess):
            nonlocal cursor
            k = cursor
            w0.append(guess)
            lbw.append(lo)
            ubw.append(hi)
            cursor += 1
            return k

        phase_pack = []  # per-phase bookkeeping for reconstruction

        for i, phase in enumerate(phases):
            sys = phase.system
            N = int(
                phase.control_intervals
                if phase.control_intervals is not None
                else control_intervals
            )
            nx, nu = sys.nx, sys.nu
            parsed = [classify_constraint(sys, c) for c in phase.constraints]

            box_initial = {("x", j): (-np.inf, np.inf) for j in range(nx)}
            box_final = {("x", j): (-np.inf, np.inf) for j in range(nx)}
            box = {("x", j): (-np.inf, np.inf) for j in range(nx)}
            box.update({("u", j): (-np.inf, np.inf) for j in range(nu)})
            t0_bounds = [-np.inf, np.inf]
            tf_bounds = [-np.inf, np.inf]
            path_constraints, initial_constraints, final_constraints = [], [], []
            for pc in parsed:
                if pc.kind == "box":
                    if pc.slot and pc.slot[0] == "t":
                        t0_bounds = [
                            max(t0_bounds[0], pc.lb),
                            min(t0_bounds[1], pc.ub),
                        ]
                        tf_bounds = [
                            max(tf_bounds[0], pc.lb),
                            min(tf_bounds[1], pc.ub),
                        ]
                    elif pc.slot:
                        lo, hi = box[pc.slot]
                        box[pc.slot] = (max(lo, pc.lb), min(hi, pc.ub))
                    else:
                        path_constraints.append(pc)
                elif pc.kind == "box_initial":
                    if pc.slot and pc.slot[0] == "t":
                        t0_bounds = [
                            max(t0_bounds[0], pc.lb),
                            min(t0_bounds[1], pc.ub),
                        ]
                    elif pc.slot:
                        lo, hi = box_initial.get(pc.slot, (-cs.inf, cs.inf))
                        box_initial[pc.slot] = (max(lo, pc.lb), min(hi, pc.ub))
                    else:
                        initial_constraints.append(pc)
                elif pc.kind == "box_final":
                    if pc.slot and pc.slot[0] == "t":
                        tf_bounds = [
                            max(tf_bounds[0], pc.lb),
                            min(tf_bounds[1], pc.ub),
                        ]
                    elif pc.slot:
                        lo, hi = box_final.get(pc.slot, (-cs.inf, cs.inf))
                        box_final[pc.slot] = (max(lo, pc.lb), min(hi, pc.ub))
                    else:
                        final_constraints.append(pc)
                elif pc.kind == "path":
                    path_constraints.append(pc)
                elif pc.kind == "initial":
                    initial_constraints.append(pc)
                elif pc.kind == "final":
                    final_constraints.append(pc)

            x_guess = [1.0] * nx
            u_guess = [1.0] * nu
            if initial_guess:
                for key, val in initial_guess.items():
                    xi = sys.state_index(key)
                    ui = sys.control_index(key)
                    if xi is not None:
                        x_guess[xi] = float(val)
                    elif ui is not None:
                        u_guess[ui] = float(val)

            pref = f"p{i}"
            T0 = cs.MX.sym(f"{pref}_T0")
            Tf = cs.MX.sym(f"{pref}_Tf")
            w += [T0, Tf]
            t0_idx.append(push_w(t0_bounds[0], t0_bounds[1], float(knots[i])))
            tf_idx.append(push_w(tf_bounds[0], tf_bounds[1], float(knots[i + 1])))
            x_idx_ph = [[] for _ in range(nx)]
            u_idx_ph = [[] for _ in range(nu)]

            g.append(Tf - T0)
            lbg.append(_MIN_PHASE_DURATION)
            ubg.append(np.inf)

            f = cs.Function(f"f_{pref}", [sys.t, sys.x, sys.u, sys.p], [sys.dx])
            has_lagrange = bool(lagrange_by_phase[i])
            lagrange_fn = None
            if has_lagrange:
                lagrange_fn = cs.Function(
                    f"L_{pref}",
                    [sys.t, sys.x, sys.u, sys.p],
                    [sum(lagrange_by_phase[i])],
                )

            def eval_expr_at(expr, t_val, x_val, u_val, _sys=sys):
                fn = cs.Function("tmp", [_sys.t, _sys.x, _sys.u, _sys.p], [expr])
                return fn(t_val, x_val, u_val, Pk)

            h_interval = (Tf - T0) / N
            Xk = cs.MX.sym(f"{pref}_X0", nx)
            w.append(Xk)
            for j in range(nx):
                lo, hi = box[("x", j)]
                lo = max(lo, box_initial[("x", j)][0])
                hi = min(hi, box_initial[("x", j)][1])
                x_idx_ph[j].append(push_w(lo, hi, x_guess[j]))

            for ic in initial_constraints:
                u_dummy0 = cs.MX.zeros(nu) if nu else cs.MX.zeros(0)
                val = eval_expr_at(ic.expr, T0, Xk, u_dummy0)
                g.append(val)
                lbg.append(ic.lb)
                ubg.append(ic.ub)

            X_all = [Xk]
            U_all = []
            Xc_all = []
            for k in range(N):
                Uk = cs.MX.sym(f"{pref}_U_{k}", nu) if nu else cs.MX.sym(f"{pref}_U_{k}", 0)
                w.append(Uk)
                for j in range(nu):
                    lo, hi = box[("u", j)]
                    u_idx_ph[j].append(push_w(lo, hi, u_guess[j]))
                U_all.append(Uk)

                Xkj = []
                for j in range(d):
                    Xkj_ = cs.MX.sym(f"{pref}_X_{k}_{j}", nx)
                    w.append(Xkj_)
                    for r in range(nx):
                        lo, hi = box[("x", r)]
                        x_idx_ph[r].append(push_w(lo, hi, x_guess[r]))
                    Xkj.append(Xkj_)

                Xk_end = D[0] * Xk
                tk = T0 + k * h_interval
                for pc in path_constraints:
                    val = eval_expr_at(pc.expr, tk, Xk, Uk)
                    g.append(val)
                    lbg.append(pc.lb)
                    ubg.append(pc.ub)
                for j in range(d):
                    xp = C[0, j + 1] * Xk
                    for r in range(d):
                        xp = xp + C[r + 1, j + 1] * Xkj[r]
                    t_cj = tk + cc.tau[j + 1] * h_interval
                    fj = f(t_cj, Xkj[j], Uk, Pk)
                    g.append(h_interval * fj - xp)
                    lbg += [0.0] * nx
                    ubg += [0.0] * nx
                    Xk_end = Xk_end + D[j + 1] * Xkj[j]
                    if has_lagrange:
                        J = J + B[j + 1] * h_interval * lagrange_fn(
                            t_cj, Xkj[j], Uk, Pk
                        )
                    for pc in path_constraints:
                        val = eval_expr_at(pc.expr, t_cj, Xkj[j], Uk)
                        g.append(val)
                        lbg.append(pc.lb)
                        ubg.append(pc.ub)

                Xk_next = cs.MX.sym(f"{pref}_X_{k+1}", nx)
                w.append(Xk_next)
                is_last = k == N - 1
                for r in range(nx):
                    lo, hi = box[("x", r)]
                    if is_last:
                        lo = max(lo, box_final[("x", r)][0])
                        hi = min(hi, box_final[("x", r)][1])
                    x_idx_ph[r].append(push_w(lo, hi, x_guess[r]))
                g.append(Xk_end - Xk_next)
                lbg += [0.0] * nx
                ubg += [0.0] * nx
                Xk = Xk_next
                X_all.append(Xk)
                Xc_all.append(Xkj)

            u_last = U_all[-1] if U_all else (
                cs.MX.zeros(nu) if nu else cs.MX.zeros(0)
            )
            for pc in path_constraints:
                val = eval_expr_at(pc.expr, Tf, Xk, u_last)
                g.append(val)
                lbg.append(pc.lb)
                ubg.append(pc.ub)
            for fc in final_constraints:
                val = eval_expr_at(fc.expr, Tf, Xk, u_last)
                g.append(val)
                lbg.append(fc.lb)
                ubg.append(fc.ub)

            for mt in mayer_by_phase[i]:
                t_val = T0 if mt.point == "0" else Tf
                x_val = X_all[0] if mt.point == "0" else X_all[-1]
                u_dummy = U_all[-1] if U_all else (
                    cs.MX.zeros(nu) if nu else cs.MX.zeros(0)
                )
                val = eval_expr_at(mt.expr, t_val, x_val, u_dummy)
                J = J + val

            phase_pack.append(
                dict(
                    sys=sys,
                    N=N,
                    d=d,
                    T0=T0,
                    Tf=Tf,
                    X_all=X_all,
                    U_all=U_all,
                    Xc_all=Xc_all,
                    nx=nx,
                    nu=nu,
                )
            )
            x_idx.append(x_idx_ph)
            u_idx.append(u_idx_ph)
            systems.append(sys)

        if self._continuity:
            for i in range(n_ph - 1):
                a, b = phase_pack[i], phase_pack[i + 1]
                g.append(a["Tf"] - b["T0"])
                lbg.append(0.0)
                ubg.append(0.0)
                nlink = min(a["nx"], b["nx"])
                if nlink == 0:
                    continue
                g.append(a["X_all"][-1][:nlink] - b["X_all"][0][:nlink])
                lbg += [0.0] * nlink
                ubg += [0.0] * nlink

        if self._sense == "max":
            J = -J

        w_sym = cs.vertcat(*w)
        g_sym = cs.vertcat(*g) if g else cs.MX.sym("g_empty", 0)
        return _CompiledMultiphase(
            nlp={"f": J, "x": w_sym, "g": g_sym},
            w_sym=w_sym,
            lbw=np.asarray(lbw, dtype=float),
            ubw=np.asarray(ubw, dtype=float),
            lbg=np.asarray(lbg, dtype=float),
            ubg=np.asarray(ubg, dtype=float),
            w0=np.asarray(w0, dtype=float),
            phase_pack=phase_pack,
            cc=cc,
            n_phases=n_ph,
            t0_idx=t0_idx,
            tf_idx=tf_idx,
            x_idx=x_idx,
            u_idx=u_idx,
            systems=systems,
        )
