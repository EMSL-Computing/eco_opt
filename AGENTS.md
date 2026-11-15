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
python -m pytest -q
python -m mkdocs build
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

## Releases and Zenodo (mandatory)

Every change pushed to `main` must be released and archived on Zenodo.
Do not leave changes on `main` without a matching Zenodo version.

1. Bump the version in `pyproject.toml`, `CITATION.cff` and
   `src/eco_opt/__init__.py` (patch for fixes, minor for features).
2. Commit, then tag `v<version>`, and push both the branch and the tag.
3. Run `python scripts/release/zenodo_release.py`. It refuses to run unless
   the tree is clean, `HEAD` is the pushed tag, and the version is new. It
   uploads `git archive` of the tag as a new version of the Zenodo record.
4. Commit the updated `scripts/release/zenodo.json` and push it. This
   bookkeeping commit is the only change exempt from the rule; it is included
   in the next release. Repository files cite only the concept DOI
   10.5281/zenodo.24239944, so no other file needs the new version DOI.
5. If a manuscript cites a specific version, update its DOI there too.

The Zenodo token stays outside the repository (`~/.config/zenodo/api.env`).

## Done

`python -m pytest -q` passes. Add tests for behavioral changes.
