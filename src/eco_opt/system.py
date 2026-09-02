"""System definitions: states, controls, and named signals.

A :class:`System` bundles a symbolic time variable, state/control/algebraic/
parameter vectors, and a user-supplied ``model`` callback that returns the
dynamics ``dx`` and a nested namespace of named signals ``y``, e.g.::

    def goddard_model(t, x, u):
        v, h, m = x[0], x[1], x[2]
        F = u[0]
        ...
        dx = cs.vertcat(dv, dh, dm)
        y = Namespace()
        y.rocket = Namespace(velocity=v, height=h, mass=m, fuelMassFlow=F)
        return dx, y

    sys = System(states=3, controls=1, model=goddard_model)
"""

from __future__ import annotations

from typing import Callable, Optional

import casadi as cs

__all__ = ["System", "Namespace"]


class Namespace:
    """A simple attribute-holding namespace, used for the ``y`` signal tree.

    Any keyword/attribute becomes an accessible field; nesting is achieved
    by assigning further ``Namespace`` instances.
    """

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __repr__(self):
        fields = ", ".join(f"{k}={v!r}" for k, v in vars(self).items())
        return f"Namespace({fields})"


class System:
    """Symbolic OCP system: states/controls/algebraics/parameters + model.

    Parameters
    ----------
    states, controls, algebraics, parameters : int
        Dimensions of each variable class (default 0 for algebraics/parameters).
    model : callable(t, x, u, z=None, p=None) -> (dx, y) or (dx, dz, y)
        User dynamics callback. Must return the state derivative ``dx`` (a
        CasADi column vector matching ``states``), optionally the algebraic
        residual ``dz``/``alg`` when ``algebraics > 0``, and a :class:`Namespace`
        tree ``y`` of named signals built from symbolic state/control/etc.
        expressions (used for constraint/objective authoring and for
        retrieving named trajectories from the solution).
    """

    def __init__(
        self,
        states: int,
        controls: int,
        model: Callable,
        algebraics: int = 0,
        parameters: int = 0,
        name: str = "sys",
    ):
        self.name = name
        self.nx = int(states)
        self.nu = int(controls)
        self.nz = int(algebraics)
        self.np = int(parameters)

        # Symbolic independent variable (time)
        self.t = cs.MX.sym(f"{name}_t")

        # Symbolic state/control/algebraic/parameter vectors
        self.x = cs.MX.sym(f"{name}_x", self.nx) if self.nx else cs.MX.sym(f"{name}_x", 0)
        self.u = cs.MX.sym(f"{name}_u", self.nu) if self.nu else cs.MX.sym(f"{name}_u", 0)
        self.z = cs.MX.sym(f"{name}_z", self.nz) if self.nz else cs.MX.sym(f"{name}_z", 0)
        self.p = cs.MX.sym(f"{name}_p", self.np) if self.np else cs.MX.sym(f"{name}_p", 0)

        self._model = model

        # Evaluate the model symbolically once to build dx, alg residual, and y.
        # `y` is optional: models may return a bare dx (or dx, dz). Default
        # to an empty Namespace so `sys.y` always exists.
        if self.nz:
            result = model(self.t, self.x, self.u, self.z, self.p)
            if isinstance(result, tuple) and len(result) == 3:
                dx, dz, y = result
            elif isinstance(result, tuple) and len(result) == 2:
                dx, dz = result
                y = Namespace()
            else:
                dx = result
                dz = None
                y = Namespace()
            self.dz = dz
        else:
            # Models declaring parameters but no algebraics still receive `p`.
            # Try (t, x, u, p) first, then (t, x, u) for parameter-free models.
            try:
                result = model(self.t, self.x, self.u, self.p) if self.np else model(self.t, self.x, self.u)
            except TypeError:
                result = model(self.t, self.x, self.u)
            if isinstance(result, tuple) and len(result) == 3:
                dx, dz, y = result
                self.dz = dz
            elif isinstance(result, tuple) and len(result) == 2:
                dx, y = result
                self.dz = None
            else:
                dx = result
                y = Namespace()
                self.dz = None
        self.dx = dx
        self.y = y

    # -- convenience accessors -------------------------------------------------
    def state_index(self, expr) -> Optional[int]:
        """Return the index r such that self.x[r] is symbolically expr, else None."""
        for i in range(self.nx):
            if cs.is_equal(self.x[i], expr, 1):
                return i
        return None

    def control_index(self, expr) -> Optional[int]:
        for i in range(self.nu):
            if cs.is_equal(self.u[i], expr, 1):
                return i
        return None

    def __repr__(self):
        return f"System(name={self.name!r}, nx={self.nx}, nu={self.nu}, nz={self.nz}, np={self.np})"
