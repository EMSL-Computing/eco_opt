# Van der Pol (`van_der_pol`)

**Question.** How does collocation behave on a smooth nonlinear
oscillator with a quadratic running cost and box bounds on the control?

## Formulation

Unforced Van der Pol is \(\ddot x-\mu(1-x^2)\dot x+x=0\). Here a
control enters as a force; \(J\) is a quadratic running cost on state
and control over a **fixed** horizon, with box bounds on \(u\).

No closed-form \(J\).

## What was done

Standard collocation instance used as a nonlinear Lagrange gate (the
role this oscillator plays in many toolboxes).

## What we found

Solver succeeds. Regression \(J^*\approx 2.888\).

![Van der Pol](../assets/examples/van_der_pol/van_der_pol.png)

## Interpretation

A smooth oscillator with a running cost is a different test from
bang–bang allocation or a free horizon. Failure here would be a
collocation bug, not an ecological finding.

## Public API

Full script: `examples/van_der_pol/van_der_pol.py`.

```python
import eco_opt
from eco_opt import Namespace


def vdp_model(t, x, u):
    x1, x2, force = x[0], x[1], u[0]
    y = Namespace(x1=x1, x2=x2, u=force)
    return [(1.0 - x2 ** 2) * x1 - x2 + force, x1], y


sys = eco_opt.System(states=2, controls=1, model=vdp_model)
ocp = eco_opt.Ocp()
ocp.min([sys.y.x1 ** 2 + sys.y.x2 ** 2 + sys.y.u ** 2])
ocp.st(
    system=sys,
    constraints=[
        (-1.0, "<=", sys.u[0], "<=", 1.0),
        (0.0, "==", eco_opt.t_0(sys.t)),
        (10.0, "==", eco_opt.t_f(sys.t)),
        (0.0, "==", eco_opt.t_0(sys.y.x1)),
        (1.0, "==", eco_opt.t_0(sys.y.x2)),
    ],
)
sol = ocp.solve(control_intervals=40, polynomial_degree=4)
```

## Run and tests

```bash
cd examples/van_der_pol
python van_der_pol.py
python plot_van_der_pol.py
python -m pytest ../../tests/test_van_der_pol.py -q
```
