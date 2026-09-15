# eco_opt

Installable Python package for posing process-based ecological ODE models
as constrained optimal-control problems (CasADi + direct collocation).

`eco_opt` is a Python implementation of the Yop modelling workflow.
See `docs/LINEAGE.md`.

## Environment

Use `python` and `pip` from the working interpreter. Do not hard-code
absolute interpreter paths in the repository.

## Package rules

- Package name is `eco_opt`.
- Must remain `pip install`-able: `import eco_opt` and `pytest` after install.
- Public API only through `eco_opt` (`System`, `Ocp`, `OcpResults`, `Phase`,
  `PhaseTrajectory`, `t_0`, `t_f`, `Namespace`, `simulate`,
  `closed_loop`, `SimulationResults`, `ClosedLoopResults`,
  `control_from_solution`). Do not teach users to import internal
  modules. Single-phase bound/parameter sweeps: `Ocp.build`,
  `Ocp.set_constraints`, `Ocp.optimize`. Algebraics:
  `System(algebraics=n)`, model returns `(dx, dz, y)`. Open-loop
  replay: `ocp.simulate(sol)` or `eco_opt.simulate(sys, x0, t, u)`.
  Closed loop (shrinking/receding horizon):
  `ocp.closed_loop(x0, dt, t_final=...)`.
- Do not hard-code example-specific solutions into the library.
- Do not cite MATLAB class names, `.m` paths, or Yop internals in library
  code, examples, or tests. Lineage belongs only in `docs/LINEAGE.md`
  and the one-line note in the package README / `__init__.py`.
- Do not add author names, machine hostnames, or absolute environment
  paths to the package.

## Tests

```bash
python -m pytest -q
```
