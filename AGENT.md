# eco_opt

Installable Python package for posing process-based ecological ODE models
as constrained optimal-control problems (CasADi + direct collocation).

`eco_opt` is a Python implementation of the Yop modelling workflow.
See `docs/LINEAGE.md`.

## Environment

Use `/home/chak803/miniforge3/envs/llm_env/bin/python` and `.../pip`.
Do not run `conda activate` in this tool environment; the interactive
hook hangs.

## Layout (umbrella)

This git repo is the **package only**. Sibling trees under
`/tahoma/emsl61845/Arjun/eco_opt_root/`:

- `eco_opt/` — this package (this git)
- `manuscript/` — MEE Applications paper (separate git)
- `reference/` — source papers and original models (not in this git)

## Package rules

- Package name is `eco_opt`.
- Must remain `pip install`-able: `import eco_opt` and `pytest` after install.
- Public API only through `eco_opt` (`System`, `Ocp`, `OcpResults`, `t_0`,
  `t_f`, `Namespace`). Do not teach users to import internal modules.
- Do not hard-code example-specific solutions into the library.
- Do not cite MATLAB class names, `.m` paths, or Yop internals in library
  code, examples, or tests. Lineage belongs only in `docs/LINEAGE.md`
  and the one-line note in the package README / `__init__.py`.

## Manuscript rules

Edit the paper in `../manuscript/`, not here. No file paths in the paper.
This paper introduces `eco_opt`: name and version it; do not cite it as a
prior publication.

## Tests

```bash
/home/chak803/miniforge3/envs/llm_env/bin/python -m pytest -q
```
