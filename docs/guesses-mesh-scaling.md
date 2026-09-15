# Guesses, mesh, scaling

Direct collocation is a nonconvex NLP. A feasible hypothesis can fail
if the starting point, the mesh, or the units are unsuitable.

## Constant initial guess

`initial_guess` maps a symbolic component to a **constant**:

```python
sol = ocp.solve(
    control_intervals=40,
    polynomial_degree=3,
    initial_guess={sys.x[0]: 1.0, sys.u[0]: 0.0},
    tf_guess=120.0,
)
```

Useful values already satisfy the boxes (non-negative pools, an
allocation in \([0,1]\)). The default fill is not a process
trajectory. A zero guess on a quantity that appears in a denominator
can make the first IPOPT iteration fail.

`initial_guess` is not a time series. `ocp.simulate(sol)` replays a
**solution**; it does not load an arbitrary simulated \(u(t)\) into
the collocation NLP. Forward integration of \(f\) with a prescribed
open-loop \(u\) (constant uptake, zero irrigation) is
`eco_opt.simulate(sys, x0, t, u)`, independent of `solve`.

## Mesh refinement

A typical first mesh is 20–40 intervals and degree 3. After
`sol.success` and a plausible \(u(t)\), `control_intervals` may be
increased (and degree, if the state is poorly resolved). First-solve
Goddard uses 59 intervals and degree 5 to match an archive; that mesh
is not required for a nine-state C–N model.

`"radau"` is often more reliable when \(u^*\) is expected to switch;
`"legendre"` matches the Goddard archive. See
[Transcription](transcription.md).

Refining the mesh **rebuilds** the NLP. The previous `sol` is not
interpolated onto the new grid. The same constant `initial_guess`
applies; `tf_guess` may be set near the coarse \(T^*\).

## Sweeps on a frozen mesh

When only **numeric** bounds or pinned parameters change, the
collocation graph can be reused:

```python
ocp.build(control_intervals=40, tf_guess=120.0, initial_guess=guess)
sol_a = ocp.optimize()
ocp.set_constraints([(new_p, "==", sys.p[0])])
sol_b = ocp.optimize()   # warm-starts from sol_a
```

`set_constraints` retunes bounds already declared in `st`. A new path
**expression** requires `st` and `build` again. After the first
successful `optimize`, the next call reuses the previous primal as
\(x_0\) (`warm_start=True`).

`build` / `optimize` applies to bound and parameter changes on a
fixed transcription, not to a change of \(f\).

## Scaling

Automatic scaling is not implemented. If pools differ by many orders
of magnitude (gC m\(^{-2}\) versus a fraction versus a thin-layer
microbial biomass), the **model** is rescaled: nondimensional groups,
or each pool divided by a typical size, so that NLP variables are
O(1).

There is no `ocp.scale` method. Path residuals in the problem
statement are likewise O(1) when a mass-balance defect is divided by
a typical flux.

Failed solves: [Troubleshooting](troubleshooting.md)
(`return_status`, print level, restoration).
