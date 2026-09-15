# API

Public names are `eco_opt.__all__`. Import them from the package
(`import eco_opt` or `from eco_opt import Namespace`). Internal
modules such as `eco_opt.system` are not the public contract.

Constraint tuples (same forms everywhere; `t_0`/`t_f` for a boundary):

```python
(lb, "<=", expr, "<=", ub)   # two-sided
(lb, "<=", expr)             # lower: lb <= expr
(ub, ">=", expr)             # upper: ub >= expr  (expr <= ub)
(eq, "==", expr)             # equality
```

`(1.0/9.0, ">=", position)` is \(x\le 1/9\). Full table:
[Constraints and time](constraints-and-time.md).

Roles of states, controls, parameters, forcings, and costs:
[Posing a problem](posing.md). Collocation and piecewise-constant
\(u\): [Transcription](transcription.md). Boxes, path constraints,
`t_0`/`t_f`, horizon: [Constraints and time](constraints-and-time.md).
`initial_guess`, `build` / `optimize`, scaling by hand:
[Guesses, mesh, scaling](guesses-mesh-scaling.md). If `solve` fails:
[Troubleshooting](troubleshooting.md).

## `System`

```python
System(states, controls, model, algebraics=0, parameters=0, name="sys")
```

Builds symbolic time `sys.t`, state `sys.x`, control `sys.u`, optional
algebraics `sys.z` and parameters `sys.p`, then calls `model` once to
get \(\dot{\mathbf{x}}\) and a `Namespace` of named signals `sys.y`.

`model` signatures:

| Dimensions | Callable | Return |
|---|---|---|
| ODE | `(t, x, u)` | `(dx, y)` or `dx` |
| ODE + parameters | `(t, x, u, p)` | `(dx, y)` or `dx` |
| DAE | `(t, x, u, z)` or `(t, x, u, z, p)` | `(dx, dz, y)` or `(dx, dz)` |

`dx` is a list or tuple of `states` expressions (`System` stacks it).
`dz` is the residual \(0=g\)
(`algebraics` rows). `y` is optional; `sys.y` is then an empty
`Namespace`. `Namespace` does not import CasADi for you: use
`import casadi as cs` when \(f\) needs `exp`, `sin`, `sign`, or
interpolants. Function list: [CasADi expressions](casadi-expressions.md).

| Attribute | Meaning |
|---|---|
| `sys.t` | time |
| `sys.x`, `sys.u` | state and control vectors |
| `sys.z`, `sys.p` | algebraics and parameters (length 0 if unused) |
| `sys.y` | named signals |
| `sys.nx`, `sys.nu`, `sys.nz`, `sys.np` | dimensions |
| `sys.name` | label |

Constraints and objectives may use `sys.x[i]` or a named field of
`sys.y`.

## `Namespace`

```python
Namespace(height=h, velocity=v, mass=m)
```

Attribute bag. Nesting is allowed (`y.rocket.height`).

## `t_0`, `t_f`

```python
t_0(expr)   # value at the initial time of the (phase) horizon
t_f(expr)   # value at the terminal time
```

Each returns a `TimepointExpression`. Use these in constraints and in
Mayer objectives. `TimepointExpression` is in `__all__` because it is
the return type; it is not constructed by user code.

## `Ocp`

```python
ocp = Ocp()
ocp.max(expr)   # or ocp.min(expr)
ocp.st(system=sys, constraints=[...])
sol = ocp.solve()
```

One of `max` or `min`. A later call replaces the earlier one. A list
of expressions is a **running** (Lagrange) cost, integrated in time.
A `t_f(expr)` term is a **terminal** (Mayer) cost. A list may mix
both.

`st(system=sys, constraints=[...])` is the single-phase form.
`st(phases=[...], continuity=True)` is the multi-phase form (below).

### `solve` and `build`

```python
sol = ocp.solve(
    control_intervals=100,
    polynomial_degree=5,
    collocation_points="legendre",
    initial_guess=None,
    ipopt_options=None,
    tf_guess=1.0,
)
```

`solve` compiles the NLP and runs IPOPT once (`warm_start=False`).
`collocation_points` is `"legendre"` or `"radau"`. `initial_guess`
maps a symbolic component (`sys.x[i]`, `sys.u[i]`, …) to a
**constant**. Extra IPOPT keys go in `ipopt_options`
([Troubleshooting](troubleshooting.md#ipopt-options)).

Bound and parameter sweeps on a frozen mesh:
[Build and optimize](build-optimize.md).

```python
ocp.build(
    control_intervals=100,
    polynomial_degree=5,
    collocation_points="legendre",
    initial_guess=None,
    tf_guess=1.0,
)
```

`build` freezes the mesh. Numeric bounds and pinned parameters may
then change without re-discretizing (single-phase):

```python
ocp.set_constraints([(new_p, "==", sys.p[0])])
sol_b = ocp.optimize(
    initial_guess=None,
    ipopt_options=None,
    tf_guess=None,
    warm_start=True,
)
```

`set_constraints` uses the same tuple syntax as `st` and updates
`ocp.constraints`. What may change after `build` (boxes, existing
path bounds, pins, IPOPT options, guess) and what requires a new
`build` is tabulated in
[Build and optimize](build-optimize.md#after-build). A **new** path
expression requires `st` and `build` again. A later `build()` keeps
the retuned bounds.

`optimize` runs IPOPT on the compiled NLP. After the first successful
call, `warm_start=True` (the default) uses the previous primal as
\(x_0\). `solve` always passes `warm_start=False`.

### Multi-phase

```python
p1 = Phase(sys_a, constraints=[...], control_intervals=20, name="ascent")
p2 = Phase(sys_b, constraints=[...], control_intervals=20, name="landing")
ocp.st(phases=[p1, p2], continuity=True)
sol = ocp.solve(tf_guess=5.0)
```

Each phase has its **own** `System` (symbols are not shared).
`control_intervals=None` inherits the value passed to `solve` /
`build`. `name` defaults to `system.name`.

`continuity=True` (default) links time and the state vector at each
junction. Algebraics and `closed_loop` are single-phase only.

### Replay

```python
sim = ocp.simulate(sol, t_grid=None, x0=None, plant=None, p=None, z0=None)
```

Open-loop integration of `plant` with `sol.u`. Defaults: `plant` is
the OCP `System`, `t_grid` is `sol.time`, `x0` is `sol.x[0]`, `p` and
`z0` are taken from `sol` when present. Same result type as
`eco_opt.simulate`.

The keyword `plant=` is the **simulated process model**
(control-theory name). It is not vegetation.

## `Phase`

```python
Phase(system, constraints=None, control_intervals=None, name=None)
```

One stage of a sequential multi-phase problem. Constraint tuples are
interpreted on that stage’s horizon \([t_0, t_f]\).

## `OcpResults`

Returned by `solve` / `optimize`.

| Attribute | Meaning |
|---|---|
| `sol.success` | IPOPT `stats["success"]` |
| `sol.return_status` | IPOPT string |
| `sol.stats` | full IPOPT stats dict |
| `sol.objective` | scalar \(J\) in the declared sense (`max` or `min`) |
| `sol.t0`, `sol.tf` | horizon |
| `sol.time` | interval boundaries, length \(N+1\) |
| `sol.x` | states on `sol.time`, shape `(N+1, nx)` |
| `sol.u` | piecewise-constant controls, shape `(N, nu)` |
| `sol.z` | algebraics on `sol.time`, shape `(N+1, nz)` |
| `sol.p` | numerical parameters |
| `sol.time_intervals` | same knots as `sol.time` (control edges) |
| `sol.time_nodes`, `sol.x_nodes` | knots **and** collocation interiors |
| `sol.N`, `sol.d` | interval count and polynomial degree |
| `sol.sys` | `System` used in the transcription |
| `sol.phases` | list of `PhaseTrajectory`; single-phase is `[self]` |
| `sol.switch_times` | junction times (empty if one phase) |

Methods:

```python
sol.signal(expr)       # expression in sys.t, x, u, z; length N+1
sol.state(i)           # sol.x[:, i]
sol.control(i)         # sol.u[:, i]
```

`sol.signal` holds the last control piece at \(T\). Plot `sol.u`
against `sol.time_intervals` as steps
([Results and plots](results-and-plots.md);
[Transcription](transcription.md)).

When every phase has the same `nx`, `sol.time` / `sol.x` / `sol.u`
are concatenated across phases (duplicate junction times dropped).
Otherwise those arrays are the **last** phase; use `sol.phases[i]`.

## `PhaseTrajectory`

One entry of `sol.phases`. Fields: `sys`, `t0`, `tf`, `time`, `x`,
`u`, `z`, `time_intervals`, `N`, `d`, `tau`, and `signal(expr)`.

## `simulate`

```python
simulate(sys, x0, t_grid, u=None, p=None, z0=None) -> SimulationResults
```

Integrates \(\dot{\mathbf{x}}=f(t,\mathbf{x},\mathbf{u})\) with CVODES
(or IDAS if there are algebraics). `t_grid` must have at least two
points. `u` may be a constant, a callable `u(t)`, or an `OcpResults`
(piecewise-constant hold of `sol.u`). Time in \(f\) is the integrator
clock.

## `closed_loop`

```python
closed_loop(
    ocp, x0, dt, t_final=None, *,
    horizon="shrinking",
    prediction_horizon=None,
    plant=None,
    p=None,
    z0=None,
    control_intervals=None,
    polynomial_degree=3,
    collocation_points="radau",
    ipopt_options=None,
    initial_guess=None,
) -> ClosedLoopResults
```

`ocp.closed_loop(x0, dt, t_final=..., **kwargs)` is the same call.

```python
sol = ocp.solve(control_intervals=14, polynomial_degree=3, tf_guess=T)
cl = ocp.closed_loop(
    x0=sol.x[0],
    dt=T / 14,
    t_final=T,
    horizon="shrinking",
)
```

`T` must be pinned in `st`. Receding: `horizon="receding"` and
`prediction_horizon=...`. A different process: `plant=` with the same
`(nx, nu)`. Posing: [Open loop vs closed loop](posing.md#open-loop-vs-closed-loop-in-a-methods-section).

At each sampling instant, the OCP is solved from the current process
state; the first control piece is applied for `dt`; the process is
integrated; the loop repeats. `horizon="shrinking"` keeps the original
terminal time. `horizon="receding"` looks `prediction_horizon` ahead
(default: the original horizon length).

`t_final` defaults to a pinned `t_f(sys.t)` in `st`. `plant` defaults
to the OCP `System` and must have the same `(nx, nu)`. Multi-phase
problems raise `NotImplementedError`.

Mesh defaults **differ from** `solve`:

| | `solve` / `build` | `closed_loop` |
|---|---|---|
| `control_intervals` | 100 | from remaining (or prediction) horizon / `dt`, unless passed |
| `polynomial_degree` | 5 | 3 |
| `collocation_points` | `"legendre"` | `"radau"` |
| IPOPT `print_level` | 0 unless `ipopt_options` overrides | 0 |
| IPOPT `max_iter` | solver default | 2000 |

Passing `control_intervals` fixes the interval count at every step.
`initial_guess` supplies extra constants; the current process state is
filled in for `sys.x[i]` when missing.

## `control_from_solution`

```python
u_fn = control_from_solution(sol)  # u_fn(t) -> 1-d ndarray, length nu
```

Piecewise-constant hold of `sol.u` on `sol.time_intervals`.

## `SimulationResults`

| Attribute | Meaning |
|---|---|
| `time` | integrator grid |
| `x` | states, `(len(time), nx)` |
| `u` | controls, `(len(time)-1, nu)` |
| `z`, `p` | algebraics on `time`; parameters |
| `t0`, `tf` | first and last times |
| `sys` | integrated `System` |

## `ClosedLoopResults`

Same time series as `SimulationResults`, plus:

| Attribute | Meaning |
|---|---|
| `horizon` | `"shrinking"` or `"receding"` |
| `statuses` | IPOPT string per step |
| `objectives` | IPOPT \(J\) of each shrinking/receding NLP (remaining horizon), not the realized integral along the applied \(u\) |
| `success` | every status is `Solve_Succeeded` or `Solved_To_Acceptable_Level` |
