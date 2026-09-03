# eco_opt

Installable Python package for posing process-based ecological ODE models
as constrained optimal-control problems (CasADi + direct collocation).

`eco_opt` is a Python implementation of the Yop modelling workflow.
See `docs/LINEAGE.md`.

## Environment

Use `python` and `pip` from the working interpreter. Do not hard-code
absolute interpreter paths in the repository.

## Layout

This git repository is the installable package. A sibling manuscript
or reference tree may exist; it is not part of this git. Do not put
filesystem paths in the paper.

## Package rules

- Package name is `eco_opt`.
- Must remain `pip install`-able: `import eco_opt` and `pytest` after install.
- Public API only through `eco_opt` (`System`, `Ocp`, `OcpResults`, `Phase`,
  `PhaseTrajectory`, `t_0`, `t_f`, `Namespace`, `simulate`,
  `closed_loop`, `SimulationResults`, `ClosedLoopResults`,
  `control_from_solution`). Do not teach users to import internal
  modules. Single-phase bound/parameter sweeps: `Ocp.build`,
  `Ocp.set_constraints`, `Ocp.optimize`. Index-1 algebraics:
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

## Manuscript rules

If a manuscript tree exists next to this package, edit it there, not
here. No file paths in the paper. A paper that introduces `eco_opt`
should name and version the software; do not cite it as a prior
publication.

## Tests

```bash
python -m pytest -q
```
