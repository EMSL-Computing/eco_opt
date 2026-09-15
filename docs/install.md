# Install

Python 3.9+, CasADi, and IPOPT (shipped with the CasADi wheel on most
platforms). MATLAB is not required at runtime. The package is not on
PyPI; install from a git clone.

```bash
git clone git@github.com:EMSL-Computing/eco_opt.git
# HTTPS: git clone https://github.com/EMSL-Computing/eco_opt.git
cd eco_opt
```

Example scripts, tests, and `work/*.mat` are in that clone. Installing
only the package directory yields `import eco_opt` without those
folders.

Use a virtual environment so `python` and `pip` are the same
interpreter. CasADi and NumPy (and, with extras, matplotlib and SciPy)
install as wheels on typical Linux, macOS, and Windows.

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install -e ".[test]"
python -c "import eco_opt; print(eco_opt.__version__, eco_opt.__all__)"
python -m pytest -q
```

`[test]` is pytest, matplotlib, and scipy. Tests import some example
plot modules at collection time; Goddard comparison plots use
`scipy.io`. Plot scripts without pytest:

```bash
python -m pip install -e ".[examples]"
```

If `import casadi` succeeds but `ocp.solve` cannot find IPOPT, install
IPOPT from the platform package manager and confirm that CasADi still
imports.

License: MIT (`LICENSE` in the repository).
