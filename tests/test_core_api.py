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
    import casadi as cs
    from eco_opt.system import Namespace

    def model(t, x, u):
        dx = cs.vertcat(-x[0] + u[0])
        y = Namespace(pos=x[0], ctrl=u[0])
        return dx, y

    sys = eco_opt.System(states=1, controls=1, model=model, name="test")
    assert sys.nx == 1
    assert sys.nu == 1
    assert sys.y.pos is not None


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
    import casadi as cs
    from eco_opt.system import Namespace

    def model(t, x, u):
        p, v = x[0], x[1]
        a = u[0]
        dx = cs.vertcat(v, a)
        y = Namespace(position=p, velocity=v, accel=a)
        return dx, y

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
