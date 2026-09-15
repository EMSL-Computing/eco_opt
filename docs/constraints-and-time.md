# Constraints and time

Every declared constraint enters the NLP. Redundant or inconsistent
inequalities can make a feasible continuous problem infeasible after
transcription.

Constraint tuples have one form:

```python
(lb, "==", expr)                 # equality
(lb, "<=", expr)                 # one-sided (lower)
(lb, "<=", expr, "<=", ub)       # two-sided
```

`op` is `"=="`, `"<="`, or `">="`. `t_0(...)` or `t_f(...)` marks a
**boundary** condition. A bare expression is a **path** condition
(the whole horizon). A one-sided infinite bound is
`(0, "<=", expr, "<=", float("inf"))`.

## Taxonomy

| Kind | Syntax | Meaning |
|---|---|---|
| Box | bound on a bare `sys.x[i]`, `sys.u[i]`, `sys.p[i]`, `sys.t` | Bound on that decision variable |
| Path | bound on a **derived** expression, no `t_0`/`t_f` | \(g(t,\mathbf{x},\mathbf{u})\le 0\) along the horizon |
| Initial / final | `t_0(expr)` / `t_f(expr)` | Point condition at \(t_0\) or \(T\) |
| Integral (isoperimetric) | auxiliary state whose derivative is the integrand; bound that state at \(T\) | A budget, e.g. \(\int u^2\,\mathrm{d}t=10\) |

A box on \(x_i\) is not the same as a path on \(x_i+u_j\). Boxes apply
to decision variables. Path expressions are residual constraints.
There is no separate integral-constraint constructor; a budget is an
auxiliary state
([iso_constraint](examples/iso_constraint.md)).

Inconsistent boxes (`mf <= m <= m0` together with `m(0)` outside that
interval) make the NLP infeasible before IPOPT starts.

## Where a path constraint is imposed

Path inequalities are applied at each **interval start** (where the
piecewise-constant \(u\) is defined), at the **collocation points**,
and at \(T\). Collocation points are interior (Legendre) or include
the right end (Radau); they omit \(\tau=0\), so the interval-start
residual is what constrains mixed \(g(\mathbf{x},\mathbf{u})\) at the
left of each piece.

A state-only box on \(\mathbf{x}(0)\) remains a `t_0` condition.

There is no `strict()` switch. Box bounds on \(\mathbf{x}\) and
\(\mathbf{u}\) already apply to every state node and every control
piece.

## Time

`eco_opt` does not treat an omitted \(t_0\) as zero or an omitted
\(T\) as a free horizon. Initial and terminal time are declared
explicitly:

```python
(0, "==", eco_opt.t_0(sys.t))              # t_0 = 0
(T, "==", eco_opt.t_f(sys.t))              # fixed horizon
(0, "<=", eco_opt.t_f(sys.t), "<=", Tmax)  # free T, upper cap
```

A **fixed** horizon is the usual seasonal plan (leaf growth over
120 d, a flowering season). A **free** horizon is a burnout or a
minimum-time problem: `t_f(sys.t)` is bounded from below by zero,
optionally from above, and `tf_guess` is the same **order of
magnitude** as \(T\). Goddard uses `tf_guess=200`; `tf_guess=1.0` on
a 40-second burn is a common cause of solver failure. See
[First solve](first-solve.md) and
[Guesses, mesh, scaling](guesses-mesh-scaling.md).

Closed loop with a shrinking window requires a remaining horizon that
still contains the planned terminal time. If the open-loop problem
had a free horizon, the closed-loop study uses a fixed \(T\). The
microbial example does that.

## Multi-phase problems

A phase is a change of **dynamics or objective**. Vegetative then
reproductive allocation with a free switch
([multiphase_cohen](examples/multiphase_cohen.md)); ascent then
landing. Continuity of time and state at the junction is
`ocp.st(phases=[p1, p2], continuity=True)`.

A bang–bang \(u(t)\) on one system remains a single-phase problem
([cohen_king_roughgarden](examples/cohen_king_roughgarden.md)).

Algebraics are not supported in multi-phase problems.
