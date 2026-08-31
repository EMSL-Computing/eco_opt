"""Timepoint-wrapped expressions: eco_opt's equivalent of Yop's t_0()/t_f().

Yop constraints can wrap a signal expression with ``t_0(expr)`` or ``t_f(expr)``
to indicate the constraint only applies at the initial or final time of the
horizon, rather than along the whole trajectory (a "path" constraint). This
module reproduces that wrapper as a lightweight Python object so ``Ocp.st``
can classify constraints exactly like ``YopOcp.parseConstraints`` /
``YopBoxConstraint`` do in MATLAB.
"""

from __future__ import annotations

import casadi as cs

__all__ = ["TimepointExpression", "t_0", "t_f"]


class TimepointExpression:
    """Wraps a CasADi expression to mark it as evaluated at t0 or tf.

    Mirrors Yop's ``t_0(expr)`` / ``t_f(expr)`` helpers used inside ``ocp.st``
    and ``ocp.min``/``ocp.max`` to select boundary values of a signal.
    """

    __slots__ = ("expr", "point")

    def __init__(self, expr, point: str):
        if point not in ("0", "f"):
            raise ValueError("point must be '0' or 'f'")
        self.expr = expr
        self.point = point

    def __repr__(self):
        return f"t_{self.point}({self.expr})"

    # Allow using TimepointExpression as a drop-in constraint operand:
    # comparisons should generally be handled explicitly by Ocp.st, but we
    # support the common patterns for convenience/debugging.
    def __le__(self, other):
        return (self, "<=", other)

    def __ge__(self, other):
        return (self, ">=", other)

    def __eq__(self, other):  # noqa: D401 - intentional operator overload
        return (self, "==", other)

    def __hash__(self):
        return id(self)


def t_0(expr):
    """Value of ``expr`` at the initial time of the horizon (Yop's t_0)."""
    return TimepointExpression(expr, "0")


def t_f(expr):
    """Value of ``expr`` at the final time of the horizon (Yop's t_f)."""
    return TimepointExpression(expr, "f")
