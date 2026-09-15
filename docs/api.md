# API

Public names are `eco_opt.__all__`. Import them from the package
(`import eco_opt` or `from eco_opt import Namespace`). Internal
modules such as `eco_opt.system` are not the public contract.

Constraint tuples have the same form everywhere: `(lb, op, expr)`
or `(lb, op, expr, op, ub)` with `op` one of `"=="`, `"<="`, `">="`.
Wrap `expr` in `t_0` or `t_f` for a boundary condition; omit the wrap
for a path constraint.

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

`model` is `(t, x, u) -> (dx, y)` or, with algebraics, `(t, x, u, z)`
(and `p` if `parameters>0`) returning `(dx, dz, y)`. `y` may be omitted.

| Attribute | Meaning |
|---|---|
| `sys.t` | time |
| `sys.x`, `sys.u` | state and control vectors |
| `sys.y` | named signals |
| `sys.nx`, `sys.nu`, `sys.nz`, `sys.np` | dimensions |

## `Namespace`

```python
Namespace(Cs=Cs, uptake=uptake, growth=G)
```

Attribute bag. Nesting is allowed (`y.rocket.height`).

## `Ocp`

```python
ocp = Ocp()
ocp.max(expr)   # or ocp.min(expr)
ocp.st(system=sys, constraints=[...])
sol = ocp.solve(control_intervals=100, polynomial_degree=5,
                collocation_points="legendre", tf_guess=1.0)
```

A list of expressions is a **running** (Lagrange) cost, integrated in
time. A `t_f(expr)` term is a **terminal** (Mayer) cost.

`collocation_points` is `"legendre"` or `"radau"`. `initial_guess` maps
a symbolic component (`sys.x[i]`, `sys.u[i]`) to a **constant**. Extra
IPOPT keys go in `ipopt_options`. Defaults: 100 intervals, degree 5,
Legendre, `tf_guess=1.0`.

### Sweeps on a frozen mesh

```python
ocp.build(control_intervals=40, tf_guess=1.0)
sol_a = ocp.optimize()
ocp.set_constraints([(new_p, "==", sys.p[0])])
sol_b = ocp.optimize()   # warm-starts from sol_a
```

`set_constraints` retunes bounds already declared in `st`. A new path
expression needs `st` + `build` again.

### Multi-phase

```python
p1 = Phase(sys_a, constraints=[...], control_intervals=20, name="veg")
p2 = Phase(sys_b, constraints=[...], control_intervals=20, name="repro")
ocp.st(phases=[p1, p2], continuity=True)
sol = ocp.solve(tf_guess=5.0)
```

`continuity=True` links time and the state vector at each junction.

### Replay and closed loop

```python
sim = ocp.simulate(sol)
cl = ocp.closed_loop(x0, dt, t_final=T, horizon="shrinking")
```

See `simulate` and `closed_loop` below. The keyword `plant=` is the
**simulated process model** (control-theory name). It is not a
botanical plant.

## `t_0`, `t_f`

```python
t_0(expr)   # value at the initial time
t_f(expr)   # value at the terminal time
```

Returns a `TimepointExpression`. Use these in constraints and in Mayer
objectives.

## `OcpResults`

Returned by `solve` / `optimize`.

| Attribute | Meaning |
|---|---|
| `sol.time` | coarse time grid (interval boundaries) |
| `sol.x` | states on that grid, shape `(N+1, nx)` |
| `sol.u` | piecewise-constant controls, shape `(N, nu)` |
| `sol.t0`, `sol.tf` | horizon |
| `sol.objective` | scalar \(J\) in the declared sense (max or min) |
| `sol.success`, `sol.return_status` | IPOPT |
| `sol.p` | numerical parameters |
| `sol.signal(expr)` | named or derived signal on `sol.time` |
| `sol.phases`, `sol.switch_times` | multi-phase only |

## `Phase`, `PhaseTrajectory`

`Phase` is one stage (its own `System`, constraints, interval count).
`PhaseTrajectory` is that stage’s slice of an `OcpResults`.

## `simulate`

```python
simulate(sys, x0, t_grid, u=None, p=None, z0=None) -> SimulationResults
```

Integrates \(\dot{\mathbf{x}}=f(t,\mathbf{x},\mathbf{u})\) with CVODES
(or IDAS if there are algebraics). `u` may be a constant, a callable
`u(t)`, or an `OcpResults`.

## `closed_loop`

```python
closed_loop(ocp, x0, dt, t_final=None, *, horizon="shrinking",
            prediction_horizon=None, plant=None, ...)
```

At each sampling instant, solve the OCP from the current state of the
simulated process, apply the first control piece for `dt`, integrate,
repeat. `horizon="shrinking"` shortens the remaining window;
`horizon="receding"` keeps a rolling prediction window.

`plant=` must have the same `(nx, nu)` as the OCP system. Default is
the OCP `System`.

## `control_from_solution`

```python
u_fn = control_from_solution(sol)  # u_fn(t) -> ndarray
```

Piecewise-constant hold of `sol.u`.

## `SimulationResults`, `ClosedLoopResults`

Time series `time`, `x`, `u` (and `z`, `p` when present). Closed loop
also stores per-step solver `statuses` and `objectives`.
