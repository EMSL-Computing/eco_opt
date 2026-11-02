# Transcription

`ocp.solve` transcribes the continuous problem into a nonlinear
program (NLP) and passes that NLP to IPOPT. It does not integrate the
ODE and then search over \(u(t)\).

The control is **piecewise constant** on `control_intervals` equal
pieces of \([t_0,T]\). The state is a polynomial of degree
`polynomial_degree` on each piece, required to satisfy
\(\dot{\mathbf{x}}=f\) at the collocation points. The running cost is
a quadrature of those pieces. Decision variables are the control on
each interval, the state at the knots and collocation nodes, the
horizon if \(T\) is free, and any parameters.

The parameterization represents \(u\) as a step on each interval. A
plot that interpolates `sol.u` as a continuous curve misrepresents
the discretization. `matplotlib.pyplot.stairs` or
`step(..., where="post")` matches the transcription. Snippets:
[Results and plots](results-and-plots.md).
`sol.time` is the coarse knot grid (interval boundaries). Interior
collocation nodes enforce the ODE; they are not `sol.time`.

`ocp.simulate(sol)` integrates the same piecewise-constant \(u^*\)
with CVODES. Agreement between the collocation state and the
integrator indicates that the mesh resolved \(f\). It is not a second
hypothesis.

## Mesh

```python
sol = ocp.solve(
    control_intervals=40,
    polynomial_degree=3,
    collocation_points="legendre",
    tf_guess=1.0,
)
```

Defaults are 100 intervals, degree 5, Legendre points, `tf_guess=1.0`.
They are **not** the Goddard archive (59 intervals, `tf_guess=200`)
or [First solve](first-solve.md) (200 intervals, Radau, `tf_guess=2`).
A typical first mesh for a new model is 20–40 intervals and degree 3.

| Setting | Role |
|---|---|
| `control_intervals` | Number of piecewise-constant \(u\) pieces |
| `polynomial_degree` | State polynomial degree on each piece (1–9) |
| `collocation_points` | `"legendre"` or `"radau"` |
| `tf_guess` | Scale of \(T\) (and of free switch times in multi-phase) |

`"legendre"` matches the Goddard archive in this package.
`"radau"` is often more reliable when \(u^*\) is discontinuous
(bang–bang controls, pulses; allocation is one example). See
[non_diff](examples/non_diff.md) and
[bang_free_tf](examples/bang_free_tf.md).

Mesh refinement can change \(J\) by a small amount. That change is
discretization error, not a change of mechanism.

## What may appear in \(f\)

The model function is a **CasADi** expression graph, not a general
Python program. Functions on symbols:
[CasADi expressions](casadi-expressions.md). A Python branch on a symbolic state or control is
invalid:

```python
# invalid — Python if on a CasADi MX
if x[0] > 0:
    rate = vmax * x[0] / (km + x[0])
```

A switch that depends on \(t\), \(\mathbf{x}\), or \(\mathbf{u}\) is
written with `casadi.if_else` (or `casadi.fmax` / `casadi.fabs`). A
Python `if` on a numeric value that is not a CasADi symbol (a Python
float, not `sys.p`) is evaluated before the graph is built.

A jump in a **forcing** (`if_else` on time) is not a jump in the
**control**. A bang–bang \(u(t)\) (allocation, thrust, a harvest
switch) is a control structure; a cutoff written in \(f\) is a
nonsmooth forcing.
See [non_diff](examples/non_diff.md).

Named diagnostics belong in the `Namespace` `y`. They are not extra
ODEs. Algebraic closures (`System(algebraics=...)`, model returns
`(dx, dz, y)`) are collocation residuals \(g=0\), not fast ODEs. See
[algebraic_qss](examples/algebraic_qss.md). Algebraics are
single-phase only.

## `sol`

| Attribute | Meaning |
|---|---|
| `sol.success`, `sol.return_status` | IPOPT |
| `sol.time` | knot times, length \(N+1\) |
| `sol.x` | states on that grid, `(N+1, nx)` |
| `sol.u` | piecewise-constant controls, `(N, nu)` |
| `sol.tf`, `sol.objective` | horizon and \(J\) in the declared sense |
| `sol.signal(expr)` | any expression in the system symbols, on `sol.time` |

`sol.objective` already includes the sign of `max` versus `min`. If
`sol.success` is false, see [Troubleshooting](troubleshooting.md).

## Further reading

Direct collocation: Betts, *Practical Methods for Optimal Control*;
Biegler, *Nonlinear Programming* (SIAM, 2010); Rao (2009). CasADi:
Andersson et al. (2019). IPOPT: Wächter and Biegler (2006). The
modeling workflow: [Yop](https://www.yoptimization.com) (Leek, 2016).
