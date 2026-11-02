# Goddard, dimensional (`goddard_benchmark`)

**Question.** What thrust history maximizes final altitude of a
vertical rocket when burnout time is free?

This is a **transcription check**, not an ecological model.
[First solve](../first-solve.md) is microbial growth.

## Formulation

| Symbol | Role | Meaning |
|---|---|---|
| \(v,h,m\) | states | velocity, altitude, mass |
| \(F\) | control | fuel mass-flow |

**Dynamics.**

\[
\dot v=\frac{F c-D(h)\,v|v|}{m}-g(h),\qquad
\dot h=v,\qquad
\dot m=-F,
\]

with exponential drag and inverse-square gravity (SI constants as in
the archived Yop run).

**Objective.** \(\max h(T)\). **Horizon.** Free.

**Constraints.** \(v(0)=h(0)=0\), \(m(0)=214.839\),
\(m(T)\ge 67.9833\), \(0\le F\le 9.525515\).

## What was done

Declared the same way as the ecological examples. Compared to
`work/goddard_benchmark/yop_goddard.mat` (59 intervals, degree 5,
Legendre).

## What we found

Altitude and \(T\) agree with the archive to relative errors on the
order of \(10^{-14}\). Control shows bang–singular–coast.

![States vs Yop](../assets/examples/goddard_benchmark/goddard_states_comparison.png)

![Control vs Yop](../assets/examples/goddard_benchmark/goddard_control_comparison.png)

## Interpretation

Before ecological trajectories are trusted, free time, bounds, and
collocation are checked on a problem with an independent reference.
Agreement here is about the transcriber, not about rockets in ecology.
This is **not** [First solve](../first-solve.md).

## Public API

Full script: `examples/goddard_benchmark/goddard_benchmark.py`.
`import casadi as cs` is for `exp`/`sign` in \(f\), not for
`Namespace`. First hour of the package is
[First solve](../first-solve.md) (microbial growth).

```python
import casadi as cs
import eco_opt
from eco_opt import Namespace


def goddard_model(t, x, u):
    v, h, m, F = x[0], x[1], x[2], u[0]
    ...
    y = Namespace()
    y.rocket = Namespace(velocity=v, height=h, mass=m, fuelMassFlow=F)
    return [dv, dh, dm], y


sys = eco_opt.System(states=3, controls=1, model=goddard_model)
ocp = eco_opt.Ocp()
ocp.max(eco_opt.t_f(sys.y.rocket.height))
ocp.st(system=sys, constraints=[...])
sol = ocp.solve(control_intervals=59, polynomial_degree=5, tf_guess=200.0)
```

## Run and tests

```bash
cd examples/goddard_benchmark
python goddard_benchmark.py
python plot_goddard_benchmark.py
python -m pytest ../../tests/test_goddard_benchmark.py -q
```

## References

Classic Goddard problem. Yop pages:
[free \(t_f\)](https://www.yoptimization.com/goddardRocketFreeTf),
[overview](https://www.yoptimization.com/goddardRocket).
The text and code here are ours.
