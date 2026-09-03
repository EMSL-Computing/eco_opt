# Install

Requires Python 3.9+, CasADi, and IPOPT (pulled in with CasADi on most
platforms). The package does not require MATLAB at runtime.

From the package root:

```bash
pip install -e ".[test]"
python -c "import eco_opt; print(eco_opt.__version__, eco_opt.__all__)"
python -m pytest -q
```

For this documentation site:

```bash
pip install -e ".[docs]"
mkdocs serve
```

If `import casadi` succeeds but `ocp.solve` cannot find IPOPT, install
IPOPT from the platform package manager and confirm CasADi still
imports.

The package is LGPL-3.0 or later.
