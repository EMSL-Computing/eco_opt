import numpy as np
import pytest

import eco_opt


def test_import_public_api():
    assert hasattr(eco_opt, "System")
    assert hasattr(eco_opt, "Ocp")
    assert hasattr(eco_opt, "OcpResults")
    assert hasattr(eco_opt, "t_0")
    assert hasattr(eco_opt, "t_f")
    assert hasattr(eco_opt, "TimepointExpression")
    assert hasattr(eco_opt, "Namespace")
    assert hasattr(eco_opt, "Phase")
    assert hasattr(eco_opt, "PhaseTrajectory")
    assert hasattr(eco_opt, "simulate")
    assert hasattr(eco_opt, "closed_loop")
    assert hasattr(eco_opt, "SimulationResults")
    assert hasattr(eco_opt, "ClosedLoopResults")
    assert hasattr(eco_opt, "control_from_solution")


def test_system_construction():
    from eco_opt import Namespace

    def model(t, x, u):
        dx = -x[0] + u[0]
        y = Namespace(pos=x[0], ctrl=u[0])
        return dx, y

    sys = eco_opt.System(states=1, controls=1, model=model, name="test")
    assert sys.nx == 1
    assert sys.nu == 1
    assert sys.y.pos is not None


def test_model_may_return_dx_as_list():
    from eco_opt import Namespace

    def model(t, x, u):
        y = Namespace(position=x[0], speed=x[1], accel=u[0])
        return [x[1], u[0]], y

    sys = eco_opt.System(states=2, controls=1, model=model, name="list_dx")
    assert sys.dx.shape == (2, 1)


def test_dx_list_length_must_match_states():
    with pytest.raises(ValueError, match="dx"):
        eco_opt.System(
            states=2,
            controls=1,
            model=lambda t, x, u: [x[0]],
            name="bad_dx",
        )


def test_timepoint_expressions():
    import casadi as cs

    x = cs.MX.sym("x")
    t0x = eco_opt.t_0(x)
    tfx = eco_opt.t_f(x)
    assert isinstance(t0x, eco_opt.TimepointExpression)
    assert t0x.point == "0"
    assert tfx.point == "f"


def test_small_solve_double_integrator():
    """Minimum-time-like small OCP: drive a double integrator to rest at x=1.

    Exercises System/Ocp/solve/OcpResults end to end on a tiny problem
    (not Goddard) to validate general applicability of the transcription.
    """
    from eco_opt import Namespace

    def model(t, x, u):
        p, v = x[0], x[1]
        a = u[0]
        y = Namespace(position=p, velocity=v, accel=a)
        return [v, a], y

    sys = eco_opt.System(states=2, controls=1, model=model, name="di")
    y = sys.y
    ocp = eco_opt.Ocp()
    # minimize integral of a^2 (minimum-energy control) -- Lagrange term
    ocp.min([y.accel ** 2])
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(y.position)),
            (0, "==", eco_opt.t_0(y.velocity)),
            (1, "==", eco_opt.t_f(y.position)),
            (0, "==", eco_opt.t_f(y.velocity)),
            (1.0, "==", eco_opt.t_0(sys.t)),
            (2.0, "==", eco_opt.t_f(sys.t)),
            (-10, "<=", y.accel, "<=", 10),
        ],
    )
    sol = ocp.solve(control_intervals=20, polynomial_degree=3)
    assert sol.success
    assert abs(sol.t0 - 1.0) < 1e-6
    assert abs(sol.tf - 2.0) < 1e-6
    assert abs(sol.x[-1, 0] - 1.0) < 1e-4
    assert abs(sol.x[-1, 1] - 0.0) < 1e-4


_INTERNAL_IMPORTS = (
    "from eco_opt.system import",
    "from eco_opt.ocp import",
    "from eco_opt.phase import",
    "from eco_opt.results import",
    "from eco_opt.simulate import",
    "from eco_opt.expressions import",
    "from eco_opt.collocation import",
    "from eco_opt.constraints import",
)


def test_example_scripts_use_public_imports():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "examples"
    offenders = []
    for path in sorted(root.rglob("*.py")):
        text = path.read_text(encoding="utf-8")
        for needle in _INTERNAL_IMPORTS:
            if needle in text:
                rel = path.relative_to(root.parent)
                offenders.append(f"{rel}: {needle}")
    assert not offenders, (
        "example scripts must import from eco_opt, not internals:\n"
        + "\n".join(offenders)
    )
