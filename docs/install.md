# Install

Requires Python 3.9+, CasADi, and IPOPT (pulled in with CasADi on most
platforms). The package does not require MATLAB at runtime.

From the package root:

```bash
pip install -e ".[test]"
python -c "import eco_opt; print(eco_opt.__version__, eco_opt.__all__)"
python -m pytest -q
```

The documentation is on Read the Docs: <https://eco-opt.readthedocs.io/>.

To preview locally, with the same interpreter you use for `eco_opt`:

```bash
python -m pip install -e ".[docs]"
python -m mkdocs serve
```

Default address is http://127.0.0.1:8000. If that port is taken:

```bash
python -m mkdocs serve -a 127.0.0.1:8001
```

Use `python -m mkdocs`, not a bare `mkdocs` from the OS package manager.
A distro `mkdocs` does not ship the Material theme and fails with
`Unrecognised theme name: 'material'`.

If `import casadi` succeeds but `ocp.solve` cannot find IPOPT, install
IPOPT from the platform package manager and confirm CasADi still
imports.

The package is LGPL-3.0 or later.
