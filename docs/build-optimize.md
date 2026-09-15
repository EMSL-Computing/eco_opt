# Build and optimize

`ocp.solve(...)` is a one-shot wrapper: it **builds** the collocation
NLP, then **optimizes** it with IPOPT. Those two steps can be run
separately. That is the way to change **numeric** bounds or pinned
parameters without assembling the CasADi graph again.

There is no separate “parameterize” call. Bounds after `build` are
changed with `set_constraints`. The initial guess is applied inside
`optimize` (`initial_guess=`, or the previous primal when
`warm_start=True`).

## After `build`

`optimize` does not write bounds. `set_constraints` writes bounds;
`optimize` then runs IPOPT.

| May change without a new `build` | How |
|---|---|
| Box on \(x_i\), \(u_i\), \(p_i\), \(z_i\) (path, or `t_0` / `t_f`) | `set_constraints`, then `optimize` |
| Numeric bounds on a path or boundary **expression already in `st`** | same |
| Pinned parameter value | same, e.g. `(3.0, "==", sys.p[0])` |
| \(t_0\) / \(T\) pins and bounds | `t_0(sys.t)` / `t_f(sys.t)` in `set_constraints` |
| IPOPT options | `optimize(ipopt_options=...)` |
| Constant initial guess | `optimize(initial_guess=...)` |
| Starting value of \(T\) (not its bounds) | `optimize(tf_guess=...)` |
| Reuse of the last primal | `optimize(warm_start=True)` (default) |

Every \(x_i\), \(u_i\), \(p_i\), \(z_i\) already has an NLP box
(default \(\pm\infty\)), so a simple bound on those symbols may be
set even if it was omitted from the original `st`. A **derived**
path residual (\(x^2\), \(x-u\), …) may be retuned only if that
expression was in `st` at `build`.

`set_constraints` is a **patch** on `ocp.constraints`: omitted
entries keep their last bounds; matched entries are replaced; a new
box on \(x,u,p,z\) is appended. A later `build()` transcribes this
updated list (mesh refine does not restore the original `st`).

| Requires `st` and/or `build` again | |
|---|---|
| Model \(f\) or `System` | `st` + `build` |
| `max` / `min` | those calls drop the compiled NLP |
| New path or boundary **expression** | `st` + `build` |
| `control_intervals`, `polynomial_degree`, `collocation_points` | `build` |
| Multi-phase `Phase` list or junction continuity | `Phase.constraints` + `build` (`set_constraints` is single-phase) |

Worked cases: path bound
[Bryson–Denham](examples/bryson_denham.md); parameter pin
[ligninolytic investment](examples/opt_lignin.md)
(`plot_figure2.py`, not the baseline `lignin_decomposition.py`).
Tests: `tests/test_nlp_reuse.py`.

## `solve`

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

| Argument | Default | Values |
|---|---|---|
| `control_intervals` | `100` | positive integer |
| `polynomial_degree` | `5` | 1–9 |
| `collocation_points` | `"legendre"` | `"legendre"`, `"radau"` |
| `initial_guess` | fill (not a process trajectory) | `dict` of symbol → constant |
| `ipopt_options` | `print_level` 0 (see [Troubleshooting](troubleshooting.md)) | IPOPT keys |
| `tf_guess` | `1.0` | scale of \(T\) and of free switch times |

`solve` always starts IPOPT cold (`warm_start=False`). Mesh meaning:
[Transcription](transcription.md). Constant guesses:
[Guesses, mesh, scaling](guesses-mesh-scaling.md).

## `build`

```python
ocp.build(
    control_intervals=100,
    polynomial_degree=5,
    collocation_points="legendre",
    initial_guess=None,
    tf_guess=1.0,
)
```

`build` transcribes \(f\), \(J\), and the constraint **expressions**
declared in `st` onto a collocation mesh. The compiled NLP is stored
on the `Ocp`. It returns `ocp` (no solution object). What may change
after this call is listed under [After `build`](#after-build).

## `set_constraints`

```python
ocp.set_constraints([
    (0.0, "==", eco_opt.t_0(sys.x[0])),
    (1.0 / 5.0, ">=", sys.y.position),
    (new_p, "==", sys.p[0]),
])
```

Same tuple syntax as `st`. Boxes on decision variables and numeric
bounds on expressions already in the NLP are written here; see
[After `build`](#after-build). A new residual requires `st` and
`build` again.

Call `build` first. Multi-phase problems cannot use `set_constraints`;
change `Phase.constraints` and `build` again.

## `optimize`

```python
sol = ocp.optimize(
    initial_guess=None,
    ipopt_options=None,
    tf_guess=None,
    warm_start=True,
)
```

Runs IPOPT on the compiled NLP. After the first successful call,
`warm_start=True` (default) uses the previous primal as \(x_0\). Pass
`initial_guess` to replace that with per-component constants.
`tf_guess` may be omitted after the first run if \(T\) is unchanged.

After `set_constraints`, pass `warm_start=False` (and `tf_guess` if
\(T\) is free). A primal that sat on the old bound can make IPOPT
exceed its iteration limit when the bound is loosened.

Returns `OcpResults`, the same type as `solve`.

## Bound sweep (frozen mesh)

Complete [Bryson–Denham](examples/bryson_denham.md) retune of the
position path bound without rebuilding:

```python
import eco_opt


def trolley_model(t, x, u):
    y = eco_opt.Namespace(position=x[0], speed=x[1], acceleration=u[0])
    return [x[1], u[0]], y


sys = eco_opt.System(states=2, controls=1, model=trolley_model, name="trolley")
y = sys.y
ocp = eco_opt.Ocp()
ocp.min([0.5 * y.acceleration ** 2])
ocp.st(
    system=sys,
    constraints=[
        (0, "==", eco_opt.t_0(y.position)),
        (1, "==", eco_opt.t_0(y.speed)),
        (1, "==", eco_opt.t_f(sys.t)),
        (0, "==", eco_opt.t_f(y.position)),
        (-1, "==", eco_opt.t_f(y.speed)),
        (1.0 / 9.0, ">=", y.position),
    ],
)

ocp.build(control_intervals=20, polynomial_degree=5, tf_guess=1.0)
sol_a = ocp.optimize(warm_start=False, tf_guess=1.0)

ocp.set_constraints([(1.0 / 5.0, ">=", y.position)])
sol_b = ocp.optimize(warm_start=False, tf_guess=1.0)

print(sol_a.success, sol_a.objective, sol_a.x[:, 0].max())
print(sol_b.success, sol_b.objective, sol_b.x[:, 0].max())
```

`sol_a` kisses \(x=1/9\) with \(J^*\approx 4\). `sol_b` uses the
looser bound \(x\le 1/5\) and a lower cost. Default
`optimize()` (`warm_start=True`) on this loosening can return
`Maximum_Iterations_Exceeded`.

The same pattern applies to a parameter pin (`sys.p[i]`) or a terminal
state bound. Changing the ODE is not this pattern: that is `st` and
`build` again.

Mesh refinement is also a new `build`. The previous `sol` is not
interpolated onto the new grid
([Guesses, mesh, scaling](guesses-mesh-scaling.md)).
