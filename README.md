# eco_opt

[![Version](https://img.shields.io/badge/version-0.1-0D7377)](https://github.com/ArjunChakrawal/eco_opt)
[![Python](https://img.shields.io/badge/python-3.9%2B-3776AB)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-2F9E44)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-Read%20the%20Docs-0D7377)](https://eco-opt.readthedocs.io/)

Process-based ecological ODEs as constrained optimal-control problems.
The right-hand side is the user’s mass balance; states, controls, and
the criterion $J$ are declared; CasADi collocation and IPOPT produce
$u^*(\cdot)$.

Python implementation of the [Yop](https://www.yoptimization.com)
modelling workflow (Leek, 2016).

```python
import eco_opt

sys = eco_opt.System(states=3, controls=1, model=goddard_model)
sig = sys.y.rocket

ocp = eco_opt.Ocp()
ocp.max(eco_opt.t_f(sig.height))
ocp.st(system=sys, constraints=[...])
sol = ocp.solve(control_intervals=59, polynomial_degree=5, tf_guess=200.0)

print(sol.success, sol.tf, sol.objective)
```

Full first script: [First solve](docs/first-solve.md)
([Goddard, dimensional](docs/examples/goddard_benchmark.md)).

## Capabilities

| Capability | |
|---|---|
| Single-phase collocation | Legendre or Radau; piecewise-constant $u$ |
| Algebraic states | `System(algebraics=n)`; model returns `(dx, dz, y)` (single-phase) |
| Multi-phase | `Ocp.st(phases=[...])`; optional continuity of time and state |
| Sweeps | `build` → `set_constraints` / `optimize` on a frozen mesh |
| Open loop | `ocp.simulate(sol)` or `eco_opt.simulate` |
| Closed loop | shrinking or receding horizon (`ocp.closed_loop`) |

Automatic scaling is not implemented. Algebraics are not supported in
multi-phase problems. The keyword `plant=` is the simulated process
model (control-theory name), not vegetation.

**Public names:** `System`, `Ocp`, `OcpResults`, `Phase`,
`PhaseTrajectory`, `t_0`, `t_f`, `Namespace`, `simulate`,
`closed_loop`, `SimulationResults`, `ClosedLoopResults`,
`control_from_solution`.

## Install

```bash
git clone https://github.com/ArjunChakrawal/eco_opt.git
cd eco_opt
python -m pip install -e ".[test]"
python -c "import eco_opt; print(eco_opt.__version__)"
python -m pytest -q
```

Python 3.9+, CasADi, and IPOPT. See [Install](docs/install.md) if the
solver is missing from the CasADi wheel.

## Documentation

<https://eco-opt.readthedocs.io/> (after the repository is public and
the RTD project is imported).

| Page | Contents |
|---|---|
| [Posing a problem](docs/posing.md) | Roles, costs, open vs closed loop |
| [Transcription](docs/transcription.md) | Piecewise-constant $u$, collocation, mesh |
| [Constraints and time](docs/constraints-and-time.md) | Boxes, path, $t_0$/$t_f$, horizon |
| [Guesses, mesh, scaling](docs/guesses-mesh-scaling.md) | Cold start, sweeps, units |
| [Examples](docs/examples/index.md) | Ecological models and textbook OCPs |
| [API](docs/api.md) | Public contract |
| [Troubleshooting](docs/troubleshooting.md) | IPOPT failures |

Preview locally:

```bash
python -m pip install -e ".[docs]"
python -m mkdocs serve
```

<http://127.0.0.1:8000> — use `python -m mkdocs` (see
[Install](docs/install.md) for the distro-MkDocs theme error and an
alternate port). `python -m mkdocs build` writes `site/`.

## Examples

| Folder | Topic |
|---|---|
| [`microbial_growth`](docs/examples/microbial_growth.md) | Finite-substrate uptake |
| [`opt_lignin`](docs/examples/opt_lignin.md) | Ligninolytic investment |
| [`root_opt`](docs/examples/root_opt.md) | Plant–soil C and N allocation |

Further examples: [docs/examples](docs/examples/index.md).

## Citation

Chakrawal, A., Schwarz, E., Mudunuru, M. K., Leek, V., and Colombi, T.
(2026). The eco_opt package: From process models to optimal strategies.
Manuscript in preparation.

Cite the software as `eco_opt` 0.1. Also cite CasADi, IPOPT, and Yop.
See [Cite](docs/cite.md) and [`CITATION.cff`](CITATION.cff).

## License

[MIT](LICENSE). Copyright (c) 2026 Arjun Chakrawal.
