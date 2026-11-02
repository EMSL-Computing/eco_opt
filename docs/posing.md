# Posing a problem

Turn a process ODE into an optimal-control problem by naming each
symbol’s **role**, then writing \(J\) and the bounds. The package does
not choose the hypothesis.

## Roles

| Role | Symbol | In `eco_opt` | Meaning |
|---|---|---|---|
| Time | \(t\) | `sys.t` | Independent variable |
| State | \(\mathbf{x}(t)\) | `sys.x` | A quantity that is integrated; requires an initial condition |
| Control | \(\mathbf{u}(t)\) | `sys.u` | A quantity the optimizer may vary in time |
| Parameter | \(\mathbf{p}\) | `sys.p` | Constant during one solve; may be retuned in a sweep |
| Algebraic | \(\mathbf{z}(t)\) | `sys.z` | Unknown on \(g=0\), not an ODE (quasi-steady closure) |
| Forcing | \(w(t)\) | inside \(f\) | Known function of time (weather, a prescribed input) |
| Diagnostic | — | `sys.y` | Named signal for \(J\) or a constraint, not extra dynamics |

`Namespace` holds diagnostics (`sys.y.rocket.height`). Nesting is
allowed. Import it from `eco_opt`. It does **not** replace CasADi
([CasADi expressions](casadi-expressions.md)):
`+`, `-`, `*`, `/` on states work without `import casadi`; `exp`,
`sin`, `sign`, and interpolants do not. A vector \(\dot{\mathbf{x}}\)
may be a list of expressions; `System` stacks it. Import `casadi as cs`
only when the right-hand side needs those functions (Goddard:
`cs.exp`, `cs.sign`). A one-state model typically needs only
`eco_opt`.

The keyword `plant=` on `closed_loop` is the **simulated process
model** (control-theory name). It is not vegetation.

## What is a control

A quantity that is a fixed number in the process model remains a
parameter or a forcing. It is \(u(t)\) only if the hypothesis allows
that quantity to vary in time in order to improve \(J\).

Examples in this package:

- irrigation;
- greenhouse heat;
- lake loading;
- harvest;
- microbial uptake;
- plant allocation;
- rocket mass-flow.

Rain in the irrigation forecast is **forcing** (\(P\equiv 0\) in the
OCP). Realised rain on the simulated water store is a process
disturbance, not a second control.

How the continuous problem becomes an NLP, and what may appear in
\(f\): [Transcription](transcription.md). Boxes, path constraints,
and whether \(T\) is fixed: [Constraints and time](constraints-and-time.md).

## Running cost vs terminal cost

Classify **per term**. An entry that is not wrapped in `t_f(...)` is
a running (Lagrange) cost, integrated from \(t_0\) to \(T\). A
`t_f(...)` term is a terminal (Mayer) cost, the value at \(T\). A
list may mix both (greenhouse: harvest at \(T\) minus heating).
`max` / `min` replace the whole criterion, so mixed \(J\) is one list.

```python
ocp.max([growth])                        # maximize ∫ G dt
ocp.min([irrigation])                    # minimize ∫ u dt
ocp.max(eco_opt.t_f(rocket.height))      # maximize h(T)
ocp.max(eco_opt.t_f(seed))               # maximize x₂(T)
ocp.min([eco_opt.t_f(-P5 * crop), P4 * heat])  # Mayer + Lagrange
```

Free horizon: bound `t_f(sys.t)` from below (and optionally above)
and pass `tf_guess` of the right order of magnitude. Time bounds are
**explicit**; omitting them does not default \(t_0=0\) or free \(T\).
See [Constraints and time](constraints-and-time.md).

## Open loop vs closed loop in a methods section

`ocp.solve` is the statement of the hypothesis: given this \(f\),
this \(J\), these bounds, and this initial state, what is
\(u^*(\cdot)\)? The **open-loop** trajectory is that statement on the
declared model. `ocp.simulate(sol)` integrates the same \(u^*\) and
checks that the collocation tape and the ODE agree.

Use **closed loop** (`ocp.closed_loop`) when the claim is about
applying the **same** hypothesis after the process has left the
forecast — rain on the water store, a wind gust, an unforecast load.
Re-planning updates the initial condition of the remaining OCP; it
does not change \(J\) or \(f\). If the simulated process is identical
to the OCP model, shrinking-horizon closed loop reprints the
open-loop tape and is a code check, not a new process result.

`horizon="shrinking"` shortens the remaining window (the usual
finite-horizon plan). `horizon="receding"` keeps a rolling prediction
window. \(T\) must be pinned (`(T, "==", t_f(sys.t))`). Signature:
[API](api.md#closed_loop).

```python
T, dt = 14.0, 1.0
sol = ocp.solve(control_intervals=14, polynomial_degree=3, tf_guess=T)
cl = ocp.closed_loop(
    x0=sol.x[0],
    dt=dt,
    t_final=T,
    horizon="shrinking",
)
# receding (rolling look-ahead):
# cl = ocp.closed_loop(
#     x0=sol.x[0], dt=dt, t_final=T,
#     horizon="receding", prediction_horizon=4.0,
# )
# unforecast process: plant=other_sys  (same nx, nu)
```

Worked scripts: [irrigation](examples/irrigation_bucket.md),
[microbial growth](examples/microbial_growth.md), and
[lignin](examples/opt_lignin.md).

Failed solves: [Guesses, mesh, scaling](guesses-mesh-scaling.md),
then [Troubleshooting](troubleshooting.md). Bound and parameter
sweeps: [Build and optimize](build-optimize.md).

## Further reading

Optimal control (statement and necessary conditions):

- Pontryagin, Boltyanskii, Gamkrelidze and Mishchenko,
  *The Mathematical Theory of Optimal Processes*.
- Bryson and Ho, *Applied Optimal Control*.
- Kirk, *Optimal Control Theory: An Introduction*.
- Liberzon,
  [*Calculus of Variations and Optimal Control Theory*](https://doi.org/10.1515/9781400842643)
  (Princeton, 2012).
- Lenhart and Workman,
  *Optimal Control Applied to Biological Models* (Chapman & Hall/CRC, 2007).

Direct collocation and the numerical stack:

- Betts,
  *Practical Methods for Optimal Control and Estimation Using
  Nonlinear Programming* (SIAM).
- Biegler,
  *Nonlinear Programming* (SIAM, 2010) — collocation of ODEs and DAEs.
- Rao (2009),
  [A survey of numerical methods for optimal control](https://doi.org/10.2514/6.2009-789).
- Andersson, Gillis, Horn, Rawlings and Diehl (2019),
  [CasADi](https://doi.org/10.1007/s11081-018-9419-0).
- Wächter and Biegler (2006),
  [IPOPT](https://doi.org/10.1007/s10107-004-0559-y).
- Leek (2016), [Yop](https://www.yoptimization.com) — modeling
  workflow that `eco_opt` follows in Python.

An ecological transcription of the same kind:

- Manzoni, Chakrawal and Ledder (2023),
  [finite-substrate microbial growth](https://doi.org/10.3389/fevo.2023.1094269).
