"""Constraint classification helpers.

A constraint tuple such as ``(0, '==', t_0(x))`` is classified based on

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

__all__ = [
    "ParsedConstraint",
    "classify_constraint",
    "is_bare_variable",
    "variable_slot",
    "merge_constraint_tuples",
]


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
    """A single normalized constraint after unpacking a constraint tuple.

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
    """Unpack a constraint tuple into (lb, mid, ub).

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
            # single-sided inequality. Exactly one side must be numeric;
            # the numeric side and operator direction determine whether it
            # is an upper or lower bound on the expression.
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


def _expr_equal(a, b) -> bool:
    try:
        return bool(cs.is_equal(a, b, 2))
    except Exception:
        return False


def _constraints_match(a: ParsedConstraint, b: ParsedConstraint) -> bool:
    if a.kind != b.kind:
        return False
    if a.kind in ("box", "box_initial", "box_final"):
        return a.slot == b.slot
    return _expr_equal(a.expr, b.expr)


def merge_constraint_tuples(system, existing, updates):
    """Patch ``existing`` with ``updates``; omitted entries keep their bounds.

    Boxes on a decision variable may be added (NLP default is ±inf). A
    derived path or boundary expression must already be in ``existing``.
    A path box on time replaces every time-slot bound so a later
    ``build()`` matches the compiled NLP.
    """
    out = list(existing)
    parsed = [classify_constraint(system, c) for c in out]
    for tup in updates:
        pc = classify_constraint(system, tup)
        if pc.kind == "box" and pc.slot and pc.slot[0] == "t":
            keep, keep_p = [], []
            for c, p in zip(out, parsed):
                if not (p.slot and p.slot[0] == "t"):
                    keep.append(c)
                    keep_p.append(p)
            out, parsed = keep, keep_p
            out.append(tup)
            parsed.append(pc)
            continue
        idx = next(
            (i for i, p in enumerate(parsed) if _constraints_match(p, pc)),
            None,
        )
        if idx is not None:
            out[idx] = tup
            parsed[idx] = pc
        elif pc.kind in ("box", "box_initial", "box_final") and pc.slot is not None:
            out.append(tup)
            parsed.append(pc)
        else:
            raise ValueError(
                "constraint is not in the compiled NLP; add it in "
                "st() and call build() again"
            )
    return out
