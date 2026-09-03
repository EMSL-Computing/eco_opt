"""eco_opt: optimal control for process-based models, built on CasADi.

This package is a Python implementation of the Yop modelling workflow:

    sys = eco_opt.System(states=3, controls=1, model=my_model)
    time = sys.t
    sig = sys.y.rocket

    ocp = eco_opt.Ocp()
    ocp.max(eco_opt.t_f(sig.height))
    ocp.st(
        system=sys,
        constraints=[
            (0, '==', eco_opt.t_0(time)),
            (0, '==', eco_opt.t_0(sig.velocity)),
            (0, '<=', sig.velocity, '<=', float('inf')),
        ],
    )
    sol = ocp.solve(control_intervals=60)

See ``Ocp`` and ``System`` for the full API.
"""

from .system import System, Namespace
from .expressions import t_0, t_f, TimepointExpression
from .ocp import Ocp
from .results import OcpResults, PhaseTrajectory
from .phase import Phase
from .simulate import (
    simulate,
    closed_loop,
    SimulationResults,
    ClosedLoopResults,
    control_from_solution,
)

__all__ = [
    "System",
    "Ocp",
    "OcpResults",
    "Phase",
    "PhaseTrajectory",
    "t_0",
    "t_f",
    "TimepointExpression",
    "Namespace",
    "simulate",
    "closed_loop",
    "SimulationResults",
    "ClosedLoopResults",
    "control_from_solution",
]

__version__ = "0.1.0"
