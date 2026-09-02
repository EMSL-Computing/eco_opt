"""eco_opt: a Python reimplementation of the MATLAB Yop optimal-control package.

Built on CasADi, ``eco_opt`` mirrors the Yop modelling workflow:

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
from .results import OcpResults

__all__ = [
    "System",
    "Ocp",
    "OcpResults",
    "t_0",
    "t_f",
    "TimepointExpression",
    "Namespace",
]

__version__ = "0.1.0"
