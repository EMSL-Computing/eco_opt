# eco_opt

Python package for posing process-based ecological ODE models as
constrained optimal-control problems. Built on CasADi and direct
orthogonal collocation; modelling grammar follows the MATLAB toolbox
[Yop](https://yoptimization.com).

This directory is the **package git**. The MEE manuscript is
`../manuscript/` (separate git). MATLAB/Yop sources are `../reference/`.

## Status

Version 0.1.0. Public API: `System`, `Ocp`, `OcpResults`, `t_0`, `t_f`,
`Namespace`. Goddard free-final-time benchmark matches the archived Yop
reference to near machine precision. Ecological examples live under
`examples/` with regression tests.

Not yet implemented: algebraic states (DAEs), multi-phase problems,
automatic scaling.

## Install

```bash
/home/chak803/miniforge3/envs/llm_env/bin/pip install -e ".[test]"
```

## Check

```bash
/home/chak803/miniforge3/envs/llm_env/bin/python -c "import eco_opt; print(eco_opt.__all__)"
/home/chak803/miniforge3/envs/llm_env/bin/python -m pytest -q
```

## Cite

Until the MEE Applications article is published, cite this software by
name and version (`eco_opt` 0.1.0). Also cite Yop, CasADi, and IPOPT.
After publication, that article is the preferred citation.

## License

GNU Lesser General Public License v3.0 or later (LGPL-3.0-or-later).
