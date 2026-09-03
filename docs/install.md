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

To preview locally:

```bash
pip install -e ".[docs]"
mkdocs serve
```

Default address is http://127.0.0.1:8000. If that port is taken:

```bash
mkdocs serve -a 127.0.0.1:8001
```

If `import casadi` succeeds but `ocp.solve` cannot find IPOPT, install
IPOPT from the platform package manager and confirm CasADi still
imports.

The package is LGPL-3.0 or later.
