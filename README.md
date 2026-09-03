# eco_opt

Python package for posing process-based ecological ODE models as
constrained optimal-control problems. Built on CasADi and direct
orthogonal collocation.

`eco_opt` is a Python implementation of the
[Yop](https://www.yoptimization.com) modelling workflow.

Version 0.1.0.

Public API: `System`, `Ocp`, `OcpResults`, `Phase`, `PhaseTrajectory`,
`t_0`, `t_f`, `Namespace`, `simulate`, `closed_loop`,
`SimulationResults`, `ClosedLoopResults`, `control_from_solution`.

Index-1 algebraic states are supported in single-phase collocation
(`System(algebraics=...)`, model returns `(dx, dz, y)`). Not yet
implemented: automatic scaling. Algebraics are not supported in
multi-phase problems. Multi-phase sequential collocation:
`Ocp.st(phases=[...])`. Bound and parameter sweeps on a fixed mesh:
`ocp.build(...)`, then `set_constraints` / `optimize` (warm start).
Open-loop replay: `ocp.simulate(sol)`. Closed loop:
`ocp.closed_loop(x0, dt, t_final=T)` (shrinking horizon by default;
`horizon='receding'` for a rolling prediction window). If the
simulated process is the OCP model, shrinking closed loop reprints
the open-loop tape. The keyword `plant=` is that simulated process
(control-theory name), not vegetation.

## Install

From the package root:

```bash
pip install -e ".[test]"
python -c "import eco_opt; print(eco_opt.__version__, eco_opt.__all__)"
python -m pytest -q
```

See the [Install](docs/install.md) page for IPOPT notes.

## Documentation

Hosted on Read the Docs: <https://eco-opt.readthedocs.io/>
(live after the repo is public and the RTD project is imported).

The docs are MkDocs Material. They do not import CasADi; a docs-only
install is enough to preview pages, math, and example write-ups.

From the package root:

```bash
pip install -e ".[docs]"
mkdocs serve
```

Then open http://127.0.0.1:8000 (Ctrl+C to stop). Edits under `docs/`
reload in the browser.

To write a static copy (same as RTD) without serving:

```bash
mkdocs build
```

Output is `site/` (gitignored). Open `site/index.html` in a browser.
`pip install -r docs/requirements.txt` is the RTD-only equivalent if
you do not want an editable package install.

Start at [First solve](docs/first-solve.md) (Goddard rocket), then
[Posing a problem](docs/posing.md) and
[Examples](docs/examples/index.md).

## Examples

Runnable modules under `examples/`. Documentation pages use the same
folder names.

**Ecological models:** `microbial_growth`, `opt_lignin`, `root_opt`,
Cohen–King allocation, Iwasa shoot–root–fruit, Lindh flowering,
multi-phase Cohen, shallow lake, greenhouse heating, irrigation
bucket.

**Traditional OCP problems:** dimensional and nondimensional Goddard,
Goddard landing, bang-bang free horizon, Bryson–Denham, isoperimetric
constraint, nondifferentiable pulse, brachistochrone, Van der Pol,
Lenhart–Workman, Wang double integrator.

**Numerical features:** quasi-steady microbes as an algebraic
closure (`algebraic_qss`). Multi-phase and open vs closed loop are
demonstrated inside the examples above.

## Cite

Until a journal article is published, cite this software by name and
version (`eco_opt` 0.1.0). Also cite CasADi, IPOPT, and Yop.

## License

GNU Lesser General Public License v3.0 or later (LGPL-3.0-or-later).
