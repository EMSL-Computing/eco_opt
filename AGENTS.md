# eco_opt

Installable Python package for posing process-based ecological ODE
models as constrained optimal-control problems (CasADi + direct
collocation). Python implementation of the Yop modelling workflow;
lineage only in `docs/LINEAGE.md` and a one-line note in the README /
`__init__.py`.

## Layout

- `src/eco_opt/` — package
- `tests/` — pytest
- `examples/` — executable scientific examples (in the git clone)
- `docs/` — MkDocs user documentation

## Environment

Use `python` and `pip` from the working interpreter. Do not hard-code
absolute interpreter paths. Not on PyPI. No conda file; no `.env`.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install -e ".[test]"
```

Plot scripts: `.[examples]`. Docs build: `pip install -r docs/requirements.txt`.

## Commands

```bash
make test                          # python -m pytest -q
make docs                          # python -m mkdocs build
```

## Public API

Import only from `eco_opt`: `System`, `Ocp`, `OcpResults`, `Phase`,
`PhaseTrajectory`, `t_0`, `t_f`, `TimepointExpression`, `Namespace`,
`simulate`, `closed_loop`, `SimulationResults`, `ClosedLoopResults`,
`control_from_solution`. Do not teach `eco_opt.system` or other
internals.

- Sweeps: `Ocp.build`, `Ocp.set_constraints`, `Ocp.optimize`
- Algebraics: `System(algebraics=n)`; model returns `(dx, dz, y)`
- Open loop: `ocp.simulate(sol)` or `eco_opt.simulate(sys, x0, t, u)`
- Closed loop: `ocp.closed_loop(x0, dt, t_final=...)`
- `dx` is a list/tuple of derivatives (or a scalar); do not use
  `cs.vertcat` in user-facing code. CasADi is for `exp` / `sin` /
  `sign` / interpolants.

## Rules

- Preserve the public API unless the user asks to change it.
- Do not hard-code example-specific solutions into the library.
- Do not cite MATLAB class names, `.m` paths, or Yop internals except
  `docs/LINEAGE.md` and the one-line README / `__init__.py` note.
- Do not add author names, machine hostnames, or absolute environment
  paths to the package.
- Examples must stay `import eco_opt` only and remain executable.
- Prefer existing abstractions before adding dependencies.
- Never commit secrets or API keys (this repo has none).
- Do not commit generated `site/`.

## Done

`make test` (or `python -m pytest -q`) passes. Add tests for
behavioral changes.
