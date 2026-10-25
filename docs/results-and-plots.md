# Results and plots

`ocp.solve` and `ocp.optimize` return an `OcpResults` object (`sol`).
Arrays are NumPy. Plotting is ordinary matplotlib; `eco_opt` does not
ship a plot method. matplotlib and scipy are the `[examples]` extra
(`python -m pip install -e ".[examples]"`; `[test]` includes them).
They are not a core dependency.

The [First solve](first-solve.md) paste ends at
`print(sol.success, sol.tf, sol.objective)`. Running
`examples/microbial_growth/microbial_growth.py` also writes the
comparison figure and prints its path. Extra panels:
`examples/microbial_growth/plot_microbial_growth.py`.
The snippets below use a three-state `sol` (Goddard) to show named
signals. Field list:
[API](api.md#ocpresults). Piecewise-constant \(u\):
[Transcription](transcription.md).

## Scalars and grids

```python
print(sol.success, sol.return_status)
print(sol.t0, sol.tf, sol.objective)
print(sol.N, sol.d)                 # intervals, polynomial degree
print(sol.time.shape, sol.x.shape, sol.u.shape)
# sol.time : (N+1,)   interval boundaries
# sol.x    : (N+1, nx) states on sol.time
# sol.u    : (N, nu)   piecewise-constant controls
```

`sol.time_intervals` is the same knot vector as `sol.time` (the edges
of the control pieces). `sol.time_nodes` / `sol.x_nodes` add the
interior collocation points.

## Named signals

```python
rocket = sol.sys.y.rocket
h = sol.signal(rocket.height)       # length N+1, on sol.time
v = sol.signal(rocket.velocity)
m = sol.signal(rocket.mass)
F = sol.control(0)                  # length N; same as sol.u[:, 0]
h_i = sol.state(1)                  # same as sol.x[:, 1]
```

`sol.signal(expr)` evaluates any CasADi expression in `sys.t`,
`sys.x`, `sys.u`, `sys.z` on `sol.time`. At \(T\) the last control
piece is held. Integer shortcuts: `sol.state(i)`, `sol.control(i)`.

Parameters: `sol.p`. Algebraics: `sol.z`, shape `(N+1, nz)`.
Multi-phase: `sol.phases[i]` is a `PhaseTrajectory` with the same
array names; `sol.switch_times` are the junctions.

## Shadow prices (costates)

Single-phase solutions also return the costate \(\lambda(t)\), the
marginal value of each state: how much \(J\) would change per extra
unit of that state at time \(t\), if the rest of the trajectory were
re-optimized. It is recovered from the NLP multipliers of the
collocation equations, so no costate equations are derived by hand.

```python
lam = sol.costate          # (N*d, nx), lambda = dJ/dx, declared sense
t_lam = sol.costate_time   # (N*d,), collocation-point times
H = sol.hamiltonian        # (N*d,), H = L + lambda' f
```

- The sign follows the declared objective. For `ocp.max(...)`, a
  positive \(\lambda_i\) means one more unit of state \(i\) raises
  \(J\); for `ocp.min(...)`, a positive \(\lambda_i\) means it raises
  the cost.
- Values sit on the collocation points, not on `sol.time`. Use
  `np.interp(sol.time, t_lam, lam[:, i])` if another grid is needed.
- `sol.hamiltonian` is a convergence check: for dynamics and running
  payoff without explicit time dependence, \(H\) is constant along an
  optimum, and \(H = 0\) when the final time is free.
- Accuracy follows the mesh. A costate that oscillates between
  collocation points usually signals an active state constraint or a
  mesh that is too coarse.
- Multi-phase results set `costate`, `costate_time`, and
  `hamiltonian` to `None`.

In the microbial-growth example the costate of substrate rises as
\(\lambda(t) = \lambda_T e^{\gamma (t - T)}\): substrate becomes
more valuable as it runs out. The optimal uptake follows from the
price, \(u = -\beta + \sqrt{e_{max}\beta(r+\beta)(\mu + 1/\lambda)}\),
so a rising \(\lambda\) means slowing uptake.

![Shadow price of substrate, the uptake it implies, and the Hamiltonian: eco_opt versus analytical](assets/examples/microbial_growth/microbial_growth_costate.png)

```bash
cd examples/microbial_growth
python plot_microbial_growth.py   # writes figs/microbial_growth_costate.png
```

`tests/test_costate.py` checks the costate against the closed form
(both point sets), a linear-quadratic and a terminal-cost problem, and
a two-state model against a finite-difference
\(\partial J^*/\partial x_0\).

## Matplotlib: states

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(3, 1, figsize=(7, 6), sharex=True)
axes[0].plot(sol.time, sol.signal(rocket.velocity))
axes[0].set_ylabel("velocity [m/s]")
axes[1].plot(sol.time, sol.signal(rocket.height))
axes[1].set_ylabel("height [m]")
axes[2].plot(sol.time, sol.signal(rocket.mass))
axes[2].set_ylabel("mass [kg]")
axes[2].set_xlabel("time [s]")
plt.show()
```

## Matplotlib: control as steps

\(u\) is constant on each interval. A line through `sol.u` against
`sol.time[:-1]` misrepresents the transcription.

```python
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.stairs(sol.u[:, 0], sol.time_intervals)
ax.set_xlabel("time [s]")
ax.set_ylabel("fuel mass flow")
plt.show()
```

`ax.stairs(values, edges)` takes `values` of length \(N\) and `edges`
of length \(N+1\). Equivalent: `ax.step(sol.time_intervals,
np.append(sol.u[:, 0], sol.u[-1, 0]), where="post")`.

## Open-loop replay

```python
sim = ocp.simulate(sol)
ax.plot(sol.time, sol.x[:, 1], label="collocation")
ax.plot(sim.time, sim.x[:, 1], "--", label="CVODES")
```

`sim.time` / `sim.x` / `sim.u` follow the same shapes as `sol` when
the default grid `sol.time` is used. Disagreement that grows with a
coarse mesh is discretisation error, not a second hypothesis.

## Saving a figure

```python
fig.savefig("goddard.png", dpi=140)
```
