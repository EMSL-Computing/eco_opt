---
name: add-eco-opt-example
description: Add or revise an executable scientific example in the eco_opt repository, including its tests and documentation. Use for repository contributions, not ordinary user modeling assistance.
---

# Add an eco_opt example

Inspect one or two existing examples with comparable mathematical features
before choosing the structure. Follow the repository-wide requirements in
`AGENTS.md`.

Create or update the parts appropriate to the requested example:

- `examples/<name>/<name>.py` for the executable model;
- `examples/<name>/README.md` for local usage and scientific context;
- `tests/test_<name>.py` for focused behavioral coverage;
- `docs/examples/<name>.md` for user documentation; and
- the documentation navigation when a new page is added.

Keep examples executable and import the package only with `import eco_opt`.
Keep model construction separate from plotting when practical. Do not encode
an example-specific solution in the library.

Test meaningful formulation and scientific invariants. Avoid brittle exact
comparisons to solver trajectories unless the solution is analytically fixed.
Keep generated figures out of tests.

## Verify the contribution

1. Run the focused example test.
2. Run `python -m pytest -q`.
3. Run `python -m mkdocs build` when documentation changes.
4. Execute the example directly when its dependencies are available.
