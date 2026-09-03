"""Sequential phases for multi-phase optimal-control problems.

A :class:`Phase` is one collocation grid with its own :class:`System`
(dynamics, controls, named signals) and its own constraint list. An
:class:`Ocp` with ``st(phases=[...])`` transcribes the phases in order
and, by default, links them by continuity of time and of the state
vector at each junction.

See ``examples/multiphase_cohen`` for a two-phase vegetative–reproductive problem.
"""

from __future__ import annotations

from typing import List, Optional, Sequence

__all__ = ["Phase"]


class Phase:
    """One stage of a sequential multi-phase OCP.

    Parameters
    ----------
    system
        Dynamics for this stage. Use a **separate** :class:`System`
        instance per phase (symbols must not be shared).
    constraints
        Constraint tuples, same form as ``Ocp.st(constraints=...)``,
        interpreted on *this* phase's horizon ``[t0, tf]``.
    control_intervals
        Collocation intervals for this phase. ``None`` inherits the
        value passed to :meth:`Ocp.solve`.
    name
        Optional label; defaults to ``system.name``.
    """

    def __init__(
        self,
        system,
        constraints: Optional[Sequence] = None,
        control_intervals: Optional[int] = None,
        name: Optional[str] = None,
    ):
        self.system = system
        self.constraints: List = list(constraints or [])
        self.control_intervals = (
            None if control_intervals is None else int(control_intervals)
        )
        self.name = name if name is not None else getattr(system, "name", "phase")

    def __repr__(self):
        n = "inherit" if self.control_intervals is None else self.control_intervals
        return (
            f"Phase(name={self.name!r}, system={self.system!r}, "
            f"n_constraints={len(self.constraints)}, control_intervals={n})"
        )
