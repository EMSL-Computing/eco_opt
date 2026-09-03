# Troubleshooting

Direct collocation is a nonlinear program. `sol.success` false, a
restoration failure, or a trajectory that ignores a bound usually
means the **NLP** is poorly scaled or poorly guessed, not that the
hypothesis is empty.

## Read the solver first

```python
print(sol.success, sol.return_status, sol.tf, sol.objective)
```

`return_status` is IPOPT’s string (`Solve_Succeeded`,
`Solved_To_Acceptable_Level`, `Restoration_Failed`,
`Maximum_Iterations_Exceeded`, …). Pass extra keys in
`ipopt_options`, for example `{"max_iter": 5000, "print_level": 5}`.

## Free horizon

If \(T\) is free, `tf_guess` must be the right **order of magnitude**
(Goddard uses `tf_guess=200`; a guess of `1.0` on a 40-second burn
often dies). Bound `t_f(sys.t)` from below by zero. An upper bound
helps if the solver runs \(T\) to infinity.

## Mesh

Start coarse (`control_intervals` 20–40), confirm a plausible
\(u(t)\), then refine. First-solve Goddard uses 59 intervals and
degree 5 to match an archive; that is not the default for a new
nine-state model.

- `"legendre"` — matches the Goddard archive in this package.
- `"radau"` — often better when \(u^*\) is discontinuous (bang–bang,
  pulses). See [non_diff](examples/non_diff.md).

`polynomial_degree` 3 is a reasonable first try; 5 is heavier.

## Bounds and guesses

Inconsistent bounds (`mf <= m <= m0` together with an initial mass
outside that interval) make the NLP infeasible before IPOPT starts.
A one-sided bound is written `(0, "<=", expr, "<=", float("inf"))`.

`initial_guess` maps a symbolic component to a constant, e.g.
`{sys.x[0]: 1.0, sys.u[0]: 0.0}`. A constant guess that already
satisfies boxes is better than the default.

Automatic **scaling** is not implemented. If states differ by many
orders of magnitude, rescale the model (work in nondimensional
groups, or divide a pool by a typical size).

## Restoration failed / iteration limit

In order, try:

1. Coarser mesh, then refine from that solution (`build` /
   `optimize` warm-start, or feed `sol` as a guess on a finer grid).
2. A better `tf_guess` and box bounds.
3. Radau if you expect switches.
4. `print_level=5` and inspect which constraint is active.
5. Whether an algebraic residual \(g=0\) is consistent (`algebraics`,
   no \(z(0)\) unless you mean it).

## Closed loop

If `closed_loop` disagrees with `solve` on the **same** process, the
sampling interval `dt` is too large or the remaining horizon is
ill-posed (free \(T\) with a shrinking window that has already
passed the planned burnout). Pin \(T\) for the closed-loop study if
the open-loop problem had a free horizon; the microbial example does
that.

If CasADi imports but the solver cannot find IPOPT, the wheel did
not ship a usable `libipopt`. Install IPOPT from the platform
package manager and confirm `import casadi` still works.
