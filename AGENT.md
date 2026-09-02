# eco_opt

Installable Python package for posing process-based ecological ODE models
as constrained optimal-control problems (CasADi + direct collocation).

Historical Yop-port specification: `docs/yop_port_spec.md`.

## Environment

Use `/home/chak803/miniforge3/envs/llm_env/bin/python` and `.../pip`.
Do not run `conda activate` in this tool environment; the interactive
hook hangs.

## Layout (umbrella)

This git repo is the **package only**. Sibling trees under
`/tahoma/emsl61845/Arjun/eco_opt_root/`:

- `eco_opt/` — this package (this git)
- `manuscript/` — MEE Applications paper (separate git)
- `reference/` — MATLAB Yop / OptLignin / literature (not in this git)

## Package rules

- Package name is `eco_opt` (not `yop`).
- Must remain `pip install`-able: `import eco_opt` and `pytest` after install.
- No MATLAB at runtime.
- Public API only through `eco_opt` (`System`, `Ocp`, `OcpResults`, `t_0`,
  `t_f`, `Namespace`). Do not teach users to import internal modules.
- Do not hard-code example-specific solutions into the library.
- Prefer general fixes over benchmark patches.

## Manuscript rules

Edit the paper in `../manuscript/`, not here. No file paths in the paper.
This paper introduces `eco_opt`: name and version it; do not cite it as a
prior publication. Cite Yop, CasADi, IPOPT, and the four source papers.

## Tests

```bash
/home/chak803/miniforge3/envs/llm_env/bin/python -m pytest -q
```
