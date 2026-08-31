"""Optimal control problem definition and direct-collocation transcription.

Mirrors MATLAB Yop's ``YopOcp`` (min/max objective, ``st()`` constraints with
t_0()/t_f() timepoint wrapping and box/path/initial/final classification) and
``YopDirectCollocation`` (Radau/Legendre orthogonal collocation on finite
elements, matching CasADi's canonical direct-collocation formulation).
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Union

import numpy as np
import casadi as cs

from .expressions import TimepointExpression
from .constraints import classify_constraint, ParsedConstraint
from .collocation import CollocationCoefficients
from .results import OcpResults

__all__ = ["Ocp"]


class Ocp:
    """Optimal control problem: objective + constraints + solve.

    Usage (mirrors Yop)::

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
    """

    def __init__(self):
        self._sense: Optional[str] = None  # 'min' or 'max'
        self._objective_terms: List = []
        self.system = None
        self.constraints: List = []

    # -- problem definition ----------------------------------------------
    def min(self, expression):
        self._sense = "min"
        self._objective_terms = list(expression) if isinstance(expression, (list, tuple)) and not isinstance(expression, TimepointExpression) else [expression]

    def max(self, expression):
        self._sense = "max"
        self._objective_terms = list(expression) if isinstance(expression, (list, tuple)) and not isinstance(expression, TimepointExpression) else [expression]

    def st(self, system=None, systems=None, constraints=None):
        self.system = system if system is not None else systems
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
        sys = self.system
        if sys is None:
            raise RuntimeError("Ocp.st(...) must be called before solve()")

        N = int(control_intervals)
        d = int(polynomial_degree)
        cc = CollocationCoefficients(d, collocation_points)
        C, D, B = cc.C, cc.D, cc.B

        nx, nu, nz, npar = sys.nx, sys.nu, sys.nz, sys.np
        if nz:
            raise NotImplementedError("algebraic states not yet supported in this transcription")

        f = cs.Function("f", [sys.t, sys.x, sys.u, sys.p], [sys.dx])

        # -- initial-guess resolution -----------------------------------------
        # Yop's default parses each state/control component to a constant
        # guess of 1.0 across the whole horizon (see YopOcp.parseInitialGuess:
        # values.State = ones(...), values.Control = ones(...)). We match that
        # default for Yop fidelity. `initial_guess` may override individual
        # components: a dict mapping a symbolic expression (sys.x[i], sys.u[i],
        # or sys.t) to a constant numeric value applied at every node of that
        # component, mirroring MATLAB's `YopInitialGuess('signals', [...],
        # 'signalValues', [...])` pattern.
        x_guess = [1.0] * nx
        u_guess = [1.0] * nu
        p_guess = [1.0] * npar
        if initial_guess:
            for key, val in initial_guess.items():
                xi = sys.state_index(key)
                ui = sys.control_index(key)
                if xi is not None:
                    x_guess[xi] = float(val)
                elif ui is not None:
                    u_guess[ui] = float(val)
                elif cs.is_equal(sys.t, key, 1):
                    # time guess is not directly used (T0/Tf come from bounds
                    # / tf_guess); accepted for MATLAB-signature parity but
                    # otherwise a no-op.
                    pass
                else:
                    raise ValueError(f"initial_guess key {key!r} is not a component of system {sys!r}")

        # Classify all user constraints up front.
        parsed = [classify_constraint(sys, c) for c in self.constraints]

        box_initial = {("x", i): (-np.inf, np.inf) for i in range(nx)}
        box_final = {("x", i): (-np.inf, np.inf) for i in range(nx)}
        box = {("x", i): (-np.inf, np.inf) for i in range(nx)}
        box.update({("u", i): (-np.inf, np.inf) for i in range(nu)})
        box.update({("p", i): (-np.inf, np.inf) for i in range(npar)})
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

        # -- Mayer / Lagrange objective terms -------------------------------
        mayer_terms = []
        lagrange_terms = []
        for term in self._objective_terms:
            if isinstance(term, TimepointExpression):
                mayer_terms.append(term)
            else:
                lagrange_terms.append(term)
        has_lagrange = len(lagrange_terms) > 0
        if has_lagrange:
            lagrange_fn = cs.Function("L", [sys.t, sys.x, sys.u, sys.p], [sum(lagrange_terms)])

        # -- NLP variable containers -----------------------------------------
        w = []
        w0 = []
        lbw = []
        ubw = []
        g = []
        lbg = []
        ubg = []

        # Parameters (shared across whole horizon)
        if npar:
            Pk = cs.MX.sym("P", npar)
            w.append(Pk)
            for i in range(npar):
                lo, hi = box[("p", i)]
                lbw.append(lo)
                ubw.append(hi)
                w0.append(p_guess[i])
        else:
            Pk = cs.MX.sym("P", 0)

        # Free initial/final time
        T0 = cs.MX.sym("T0")
        Tf = cs.MX.sym("Tf")
        w += [T0, Tf]
        lbw += [t0_bounds[0], tf_bounds[0]]
        ubw += [t0_bounds[1], tf_bounds[1]]
        w0 += [0.0, tf_guess]

        h_interval = (Tf - T0) / N

        # Initial state
        Xk = cs.MX.sym("X0", nx)
        w.append(Xk)
        for i in range(nx):
            lo, hi = box[("x", i)]
            lo = max(lo, box_initial[("x", i)][0])
            hi = min(hi, box_initial[("x", i)][1])
            lbw.append(lo)
            ubw.append(hi)
            w0.append(x_guess[i])

        # Initial-time equality/path constraints on derived expressions
        t0_expr_fn_cache = {}

        def eval_expr_at(expr, t_val, x_val, u_val):
            fn = cs.Function("tmp", [sys.t, sys.x, sys.u, sys.p], [expr])
            return fn(t_val, x_val, u_val, Pk)

        for ic in initial_constraints:
            u_dummy0 = cs.MX.zeros(nu) if nu else cs.MX.zeros(0)
            val = eval_expr_at(ic.expr, T0, Xk, u_dummy0)
            # initial constraints generally shouldn't depend on u; if they do this is approximate
            g.append(val)
            lbg.append(ic.lb)
            ubg.append(ic.ub)

        X_all = [Xk]  # state at start of each interval (for trajectory reconstruction)
        U_all = []
        Xc_all = []  # collocation-point states per interval
        t_grid_starts = []

        J = 0.0
        t_running = T0

        for k in range(N):
            Uk = cs.MX.sym(f"U_{k}", nu) if nu else cs.MX.sym(f"U_{k}", 0)
            w.append(Uk)
            for i in range(nu):
                lo, hi = box[("u", i)]
                lbw.append(lo)
                ubw.append(hi)
                w0.append(u_guess[i])
            U_all.append(Uk)

            Xkj = []
            for j in range(d):
                Xkj_ = cs.MX.sym(f"X_{k}_{j}", nx)
                w.append(Xkj_)
                for i in range(nx):
                    lo, hi = box[("x", i)]
                    lbw.append(lo)
                    ubw.append(hi)
                    w0.append(x_guess[i])
                Xkj.append(Xkj_)

            Xk_end = D[0] * Xk
            tk = T0 + k * h_interval
            for j in range(d):
                xp = C[0, j + 1] * Xk
                for r in range(d):
                    xp = xp + C[r + 1, j + 1] * Xkj[r]

                tau_j = cc.tau[j + 1]
                t_cj = tk + tau_j * h_interval

                fj = f(t_cj, Xkj[j], Uk, Pk)
                g.append(h_interval * fj - xp)
                lbg += [0.0] * nx
                ubg += [0.0] * nx

                Xk_end = Xk_end + D[j + 1] * Xkj[j]

                if has_lagrange:
                    J = J + B[j + 1] * h_interval * lagrange_fn(t_cj, Xkj[j], Uk, Pk)

                # path constraints evaluated at every collocation node
                for pc in path_constraints:
                    val = eval_expr_at(pc.expr, t_cj, Xkj[j], Uk)
                    g.append(val)
                    lbg.append(pc.lb)
                    ubg.append(pc.ub)

            Xk_next = cs.MX.sym(f"X_{k+1}", nx)
            w.append(Xk_next)
            is_last = k == N - 1
            for i in range(nx):
                lo, hi = box[("x", i)]
                if is_last:
                    lo = max(lo, box_final[("x", i)][0])
                    hi = min(hi, box_final[("x", i)][1])
                lbw.append(lo)
                ubw.append(hi)
                w0.append(x_guess[i])

            g.append(Xk_end - Xk_next)
            lbg += [0.0] * nx
            ubg += [0.0] * nx

            Xk = Xk_next
            X_all.append(Xk)
            Xc_all.append(Xkj)

        # Final-time equality/path constraints on derived expressions
        for fc in final_constraints:
            u_last = U_all[-1] if U_all else (cs.MX.zeros(nu) if nu else cs.MX.zeros(0))
            val = eval_expr_at(fc.expr, Tf, Xk, u_last)
            g.append(val)
            lbg.append(fc.lb)
            ubg.append(fc.ub)

        # Mayer objective term(s): evaluated at (Tf, Xf) or (T0, X0)
        for mt in mayer_terms:
            t_val = T0 if mt.point == "0" else Tf
            x_val = X_all[0] if mt.point == "0" else X_all[-1]
            u_dummy = U_all[-1] if U_all else (cs.MX.zeros(nu) if nu else cs.MX.zeros(0))
            val = eval_expr_at(mt.expr, t_val, x_val, u_dummy)
            J = J + val

        if self._sense == "max":
            J = -J

        w = cs.vertcat(*w)
        g = cs.vertcat(*g) if g else cs.MX.sym("g_empty", 0)

        nlp = {"f": J, "x": w, "g": g}
        opts = {"ipopt": dict(ipopt_options or {})}
        opts["ipopt"].setdefault("print_level", 0)
        opts.setdefault("print_time", 0)
        solver = cs.nlpsol("solver", "ipopt", nlp, opts)

        sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
        stats = solver.stats()

        return OcpResults(
            sys=sys,
            sol=sol,
            stats=stats,
            N=N,
            d=d,
            npar=npar,
            tau=cc.tau,
            T0_sym=T0,
            Tf_sym=Tf,
            Pk_sym=Pk,
            X_all=X_all,
            U_all=U_all,
            Xc_all=Xc_all,
            w=w,
            objective_sense=self._sense,
        )
