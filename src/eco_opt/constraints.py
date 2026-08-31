"""Constraint classification helpers.

Mirrors Yop's ``YopBoxConstraint`` / ``YopPathConstraint`` /
``YopInitialConstraint`` / ``YopFinalConstraint`` predicates
(``Constraints/*.m``): a constraint tuple such as
``(0, '==', t_0(x))`` is classified based on

* whether the middle operand is wrapped in ``t_0(...)``/``t_f(...)``
  (boundary) or bare (applies along the whole horizon => "path"), and
* whether the wrapped/bare operand is itself a plain system variable
  (state/control/algebraic/parameter/time component) => "box" constraint
  (simple bound), or a derived expression => general nonlinear constraint.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import casadi as cs

from .expressions import TimepointExpression

__all__ = ["ParsedConstraint", "classify_constraint", "is_bare_variable", "variable_slot"]


def _all_components(system):
    """Yield (kind, index, symbol) for every scalar component of a system."""
    yield ("t", None, system.t)
    for i in range(system.nx):
        yield ("x", i, system.x[i])
    for i in range(system.nu):
        yield ("u", i, system.u[i])
    for i in range(system.nz):
        yield ("z", i, system.z[i])
    for i in range(system.np):
        yield ("p", i, system.p[i])


def variable_slot(system, expr) -> Optional[tuple]:
    """If expr is symbolically identical to a bare system variable component,
    return (kind, index); otherwise None."""
    for kind, idx, sym in _all_components(system):
        try:
            if cs.is_equal(sym, expr, 5):
                return (kind, idx)
        except Exception:
            continue
    return None


def is_bare_variable(system, expr) -> bool:
    return variable_slot(system, expr) is not None


@dataclass
class ParsedConstraint:
    """A single normalized constraint after unpacking a Yop-style tuple.

    kind: one of 'box', 'box_initial', 'box_final', 'path', 'initial', 'final'
    lb, ub: numeric bounds (float, possibly +/-inf)
    expr: the (unwrapped) CasADi expression the bound applies to
    slot: (kind, index) if expr is a bare variable component, else None
    """

    kind: str
    lb: float
    ub: float
    expr: Any
    slot: Optional[tuple]


def _unpack_tuple(t):
    """Unpack a Yop-style constraint tuple into (lb, mid, ub).

    Supported forms:
        (lb, '<=', mid, '<=', ub)
        (lb, '==', mid)
        (mid, '==', ub)  (rare, symmetric equality)
    """
    if len(t) == 5:
        lb, op1, mid, op2, ub = t
        if op1 not in ("<=", ">=") or op2 not in ("<=", ">="):
            raise ValueError(f"unsupported bound operators in {t}")
        # normalize to lb <= mid <= ub
        if op1 == ">=":
            lb, mid = mid, lb
        if op2 == ">=":
            mid, ub = ub, mid
        return float(lb), mid, float(ub)
    elif len(t) == 3:
        a, op, b = t
        if op == "==":
            # one side is numeric, the other is the expression
            if isinstance(a, (int, float)):
                return float(a), b, float(a)
            elif isinstance(b, (int, float)):
                return float(b), a, float(b)
            else:
                raise ValueError(f"cannot classify equality constraint {t}: need one numeric side")
        elif op in ("<=", ">="):
            # single-sided inequality: e.g. {1/9 '>=' trolley.position} means
            # trolley.position <= 1/9 (upper bound only, lb stays -inf), or
            # {expr '<=' bound} meaning expr <= bound. Exactly one side must
            # be numeric; the numeric side + operator direction determines
            # whether it's an upper or lower bound on the expression.
            if isinstance(a, (int, float)) and not isinstance(b, (int, float)):
                num, expr = float(a), b
                # a op b  =>  num op expr
                if op == ">=":
                    # num >= expr  =>  expr <= num  (upper bound)
                    return -float("inf"), expr, num
                else:
                    # num <= expr  =>  expr >= num  (lower bound)
                    return num, expr, float("inf")
            elif isinstance(b, (int, float)) and not isinstance(a, (int, float)):
                expr, num = a, float(b)
                # expr op num
                if op == "<=":
                    # expr <= num  (upper bound)
                    return -float("inf"), expr, num
                else:
                    # expr >= num  (lower bound)
                    return num, expr, float("inf")
            else:
                raise ValueError(f"cannot classify single-sided constraint {t}: need exactly one numeric side")
        else:
            raise ValueError(f"unsupported 3-tuple operator in {t}")
    else:
        raise ValueError(f"unsupported constraint tuple length {len(t)}: {t}")


def classify_constraint(system, constraint_tuple) -> ParsedConstraint:
    lb, mid, ub = _unpack_tuple(constraint_tuple)

    if isinstance(mid, TimepointExpression):
        point = mid.point
        inner = mid.expr
        slot = variable_slot(system, inner)
        if slot is not None:
            kind = "box_initial" if point == "0" else "box_final"
        else:
            kind = "initial" if point == "0" else "final"
        return ParsedConstraint(kind=kind, lb=lb, ub=ub, expr=inner, slot=slot)
    else:
        slot = variable_slot(system, mid)
        kind = "box" if slot is not None else "path"
        return ParsedConstraint(kind=kind, lb=lb, ub=ub, expr=mid, slot=slot)
