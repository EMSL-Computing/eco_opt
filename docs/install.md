# Install

Requires Python 3.9+, CasADi, and IPOPT (pulled in with CasADi on most
platforms). The package does not require MATLAB at runtime.

A wheel or `python -m pip install eco_opt` installs the library only
(`import eco_opt`). Example scripts, tests, and `work/*.mat` live in
the **git clone**, not in the wheel.

A virtual environment is cheap (`python -m venv` is seconds). The
wait is downloading wheels (CasADi, NumPy; extras add matplotlib and
SciPy). Those wheels do not compile on typical Linux, macOS, or
Windows. Use a venv so `python` and `pip` are the same interpreter
(a distro `pip` often is not).

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install -e ".[test]"
python -c "import eco_opt; print(eco_opt.__version__, eco_opt.__all__)"
python -m pytest -q
```

`[test]` is pytest plus matplotlib and scipy (tests import some
example plot modules at collection time; Goddard comparison plots
use `scipy.io`). Plot scripts alone:

```bash
python -m pip install -e ".[examples]"
```

The documentation is on Read the Docs: <https://eco-opt.readthedocs.io/>.

To preview locally, with the same interpreter as `eco_opt`:

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

The package is MIT licensed. See the `LICENSE` file in the repository.
