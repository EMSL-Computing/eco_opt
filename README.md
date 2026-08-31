# eco_opt

A CasADi-based Python reimplementation of the MATLAB [Yop](https://yoptimization.com)
optimal-control toolbox. `eco_opt` preserves Yop's modelling workflow
(`YopSystem` / `YopOcp` / `t_0()` / `t_f()` / direct collocation via
`nlpsol`+IPOPT) while being an idiomatic, installable Python package.

## Status

* Core package implemented: `System`, `Ocp`, `OcpResults`, `t_0`, `t_f`,
  `TimepointExpression`.
* Transcription: general Radau/Legendre orthogonal collocation (parametrized
  by `control_intervals` and `polynomial_degree`, matching Yop's
  `YopCollocationCoefficients` / `YopDirectCollocation`), free or fixed
  initial/final time, box/path/initial/final constraint classification
  mirroring Yop's `Constraints/*.m` predicates.
* **Mandatory first gate passed**: the Goddard rocket free-final-time
  benchmark reproduces the MATLAB Yop reference trajectory
  (`work/goddard_benchmark/yop_goddard.mat`) to near machine precision.
* Not yet implemented (follow-up work): algebraic states (DAEs),
  multi-phase/multi-system OCPs, scaling, bangFreeTf / bryson / fixed-tf
  Goddard / isoConstraint / nonDiff examples.

## Install

```bash
/home/chak803/miniforge3/envs/llm_env/bin/pip install -e ".[test]"
```

## Quick check

```bash
/home/chak803/miniforge3/envs/llm_env/bin/python -c "import eco_opt; print(eco_opt.__all__)"
```

## Public API

```python
import eco_opt
from eco_opt.system import Namespace

def model(t, x, u):
    # x, u are CasADi MX vectors; return (dx, y) where y is a Namespace tree
    # of named signals (mirrors Yop's y.rocket.velocity struct access).
    dx = ...
    y = Namespace()
    y.rocket = Namespace(velocity=..., height=..., mass=..., fuelMassFlow=...)
    return dx, y

sys = eco_opt.System(states=3, controls=1, model=model)
time = sys.t
rocket = sys.y.rocket

ocp = eco_opt.Ocp()
ocp.max(eco_opt.t_f(rocket.height))          # Mayer term (t_0/t_f wrapped)
# ocp.min([running_cost_expr])               # Lagrange term (bare expr, list)
ocp.st(
    system=sys,
    constraints=[
        (0, "==", eco_opt.t_0(time)),          # initial-time box constraint
        (0, "==", eco_opt.t_0(rocket.velocity)),
        (0, "<=", eco_opt.t_f(time), "<=", float("inf")),  # free tf, tf>=0
        (0, "<=", rocket.velocity, "<=", float("inf")),    # path box constraint
        (0, "<=", rocket.fuelMassFlow, "<=", 9.525515),    # control bound
    ],
)
sol = ocp.solve(control_intervals=59, polynomial_degree=5,
                collocation_points="legendre")

print(sol.success, sol.return_status, sol.objective, sol.tf)
v = sol.signal(rocket.velocity)   # sampled on sol.time (coarse interval grid)
h = sol.state(1)                  # by state index
F = sol.control(0)                # piecewise-constant control trajectory
```

### Constraint classification (mirrors `Constraints/*.m`)

A constraint tuple `(lb, '<=', expr, '<=', ub)` or `(val, '==', expr)` is
classified by `eco_opt.constraints.classify_constraint`:

* `expr` wrapped in `t_0(...)`/`t_f(...)` **and** a bare system variable →
  `box_initial` / `box_final` (simple bound at one endpoint, e.g. `m(0)=m0`).
* `expr` wrapped in `t_0(...)`/`t_f(...)` but a *derived* expression → general
  `initial` / `final` nonlinear constraint.
* bare `expr` that is a plain system variable → `box` (applies at every
  collocation node as a simple variable bound — matches Yop's box handling).
* bare `expr` that is a derived expression → `path` constraint, enforced at
  every collocation node.

### Results / plotting equivalent

Yop's `sol.plot(time, signal)` / `sol.stairs(time, signal)` are replaced by
Pythonic sampling: `OcpResults.signal(expr)` evaluates any CasADi expression
built from `sys.x`/`sys.u` on the solution's coarse (interval-boundary) time
grid and returns a NumPy array; combine with `matplotlib.pyplot.plot(sol.time, v)`
yourself for plotting.

## Examples

```bash
/home/chak803/miniforge3/envs/llm_env/bin/python examples/goddard_benchmark.py
```

Reproduces the Goddard rocket maximum-ascent, free-final-time benchmark
(`yoptimization/pages/benchmark/goddard_benchmark.md`), 59 control intervals,
degree-5 Legendre collocation — matching Yop's exact configuration.

### Ecology/publication examples

* `examples/yop_microbial_growth/` — microbial growth optimization for a
  finite food source (Manzoni et al. 2022), free terminal time, port of
  `examples/OptLignin/YOP benchmark/Manzoni et al. 2022/`. See its own
  README for problem statement, run instructions, and validation.
* `examples/opt_lignin_figure2/` — optimal ligninolytic activity during
  litter decomposition (Chakrawal et al. 2024, OptLignin Figure 2), fixed
  terminal time, parameterized model, port of
  `examples/OptLignin/OptLignin1.0/opt_con.m` + `Figure2.m`. See its own
  README for problem statement, run instructions, and validation.

## Tests

```bash
/home/chak803/miniforge3/envs/llm_env/bin/python -m pytest -q
```

Covers: package import & public API surface, `System`/`Ocp` construction,
`TimepointExpression` semantics, a small double-integrator solve/result
retrieval, and the Goddard benchmark regression test (objective, final time,
and interpolated state/control trajectories numerically compared against
`work/goddard_benchmark/yop_goddard.mat` with documented tolerances — see
`tests/test_goddard_benchmark.py` docstring for the exact tolerances and
rationale).

## Numerical validation summary (Goddard benchmark)

Computed with `control_intervals=59`, `polynomial_degree=5`, `'legendre'`
points (same as the MATLAB Yop benchmark configuration):

| quantity            | eco_opt         | yop_goddard.mat | abs/rel diff        |
|----------------------|-----------------|-----------------|----------------------|
| objective (h(tf))    | 161413.472864   | 161413.472864   | 1.16e-9 / 7.2e-15    |
| final time tf        | 212.954467      | 212.954467      | 6.0e-13 / 2.8e-15    |
| velocity RMS (interp) | —              | —               | 6.5e-12              |
| height RMS (interp)   | —              | —               | 5.6e-10              |
| mass RMS (interp)     | —              | —               | 6.0e-7               |
| control RMS (interp)  | —              | —               | 1.2e-8               |

Agreement is near machine precision because eco_opt reproduces the identical
Radau/Legendre collocation scheme, node count, and constraint formulation
used by Yop's `YopDirectCollocation`/`YopCollocationCoefficients`.

## Known limitations

* No algebraic-state (DAE) support yet (`System(algebraics=...)` raises
  `NotImplementedError` inside `Ocp.solve`).
* Single-system/single-phase OCPs only (no multi-phase stitching yet, though
  `YopOcp.optimize` in MATLAB supports it).
* No automatic scaling (`YopOcp.Scaling`) — not needed for Goddard, but will
  matter for stiffer problems later.
* `initial`/`final` nonlinear (non-box) constraints assume they don't
  meaningfully depend on the control (a zero-vector dummy control is
  substituted); fine for Goddard (no such constraints used) but should be
  revisited before those constraint kinds are exercised for real.
