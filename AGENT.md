# Objective

Reimplement the MATLAB Yop optimal-control package as a functioning, installable Python package named `eco_opt`.

Do not produce a `.zip` archive of scripts as the deliverable. The deliverable must be a usable Python package that can be installed, imported, tested, and used by downstream Python code.

Primary requirement: **behavioral and numerical compatibility with the original MATLAB Yop implementation**, not merely similar-looking Python code.

Keep responses concise. Spend tokens on implementation, testing, debugging, and package quality rather than explanation.

# Required deliverable

Create a proper Python package with:

* package name: `eco_opt`
* importable modules and public APIs
* `pyproject.toml` or equivalent packaging configuration
* installable package using `pip install`
* documented dependencies
* automated tests
* examples or benchmark entry points
* no requirement for users to unpack files manually
* no dependence on MATLAB at runtime

The package must support usage such as:

```python
import eco_opt
```

and expose the implemented Yop-compatible functionality through stable Python modules and documented APIs.

Do not rename the package to `yop`. The Python package must be called `eco_opt`, while preserving Yop-compatible concepts and behavior.

# Environment

Always run commands inside:

```bash
conda activate /home/chak803/miniforge3/envs/llm_env/
```

Use subagents for package installation or other independent auxiliary tasks when useful.

If context or token limits become restrictive, delegate self-contained implementation or debugging tasks to Codex GPT Luna rather than stopping unfinished.

# Sources of truth

Use these local resources before making assumptions:

* `Arjun/eco_opt/yoptimalcontrol.com/` — offline documentation for the MATLAB Yop package.
* Existing MATLAB source code — authoritative reference for API semantics and algorithms.
* `Arjun/eco_opt/yoptimization/pages/` — examples and benchmarks that define expected behavior.

Do not redesign the mathematical formulation merely because another Python formulation is easier.

When behavior is unclear, inspect the MATLAB implementation and documentation first.

# Package implementation goal

Create an idiomatic, maintainable Python package named `eco_opt` exposing functionality equivalent to the MATLAB Yop package.

Preserve, where practical:

* optimal-control problem semantics
* variable, state, and control definitions
* bounds and constraints
* derivative handling
* objective and integral definitions
* time-domain behavior
* transcription and discretization semantics
* solver behavior and tolerances
* initial guesses
* returned trajectories and objective values
* diagnostic and convergence information

Python syntax and module organization may differ where necessary, but mathematical meaning and numerical behavior must remain equivalent.

Avoid hard-coding benchmark-specific solutions.

# Public API requirements

Design and document a coherent Python API for:

* defining optimal-control problems
* declaring states, controls, parameters, and time variables
* defining dynamics
* defining running and terminal costs
* defining boundary, path, and terminal constraints
* specifying bounds and initial guesses
* selecting discretization and solver options
* solving the problem
* retrieving optimized trajectories, controls, parameters, objective values, and solver status

Use Python-native conventions, but keep the API conceptually close enough to Yop that MATLAB examples can be translated systematically.

Expose public symbols deliberately through `eco_opt.__init__` or documented submodules. Avoid requiring users to import internal implementation modules.

# Packaging requirements

The repository must contain a standard Python package layout, for example:

```text
eco_opt/
├── pyproject.toml
├── README.md
├── src/
│   └── eco_opt/
│       ├── __init__.py
│       ├── ...
│       └── tests/
```

The exact layout may differ, but the following must work from a clean environment:

```bash
pip install .
python -c "import eco_opt"
pytest
```

Declare runtime and test dependencies explicitly.

Do not rely on files being present in the current working directory after installation.

# Development strategy

Work incrementally.

1. Inspect and map the MATLAB package architecture.
2. Inspect the documentation and benchmark examples.
3. Define the `eco_opt` package structure and public API.
4. Identify the minimum functionality required by the Goddard benchmark.
5. Implement that functionality as reusable package components.
6. Install the package locally with `pip install .`.
7. Run the benchmark through the installed package.
8. Diagnose discrepancies against MATLAB/reference results.
9. Iterate until the benchmark passes.
10. Expand functionality for additional examples.
11. Add regression tests for every validated feature.

Do not create a collection of one-off scripts. Every implemented capability must belong to the reusable `eco_opt` package.

Do not attempt to port the entire package before establishing a working vertical slice, but ensure that the vertical slice is implemented as production package code rather than benchmark-only code.

# Mandatory first gate — Goddard benchmark

Start exclusively with:

`Arjun/eco_opt/yoptimization/pages/benchmark/goddard_benchmark.zip`

Reproduce the benchmark numerically using the installed `eco_opt` package.

The benchmark must be runnable through a Python module, test, or documented example, rather than only through an ad hoc script.

Compare important outputs, including where available:

* objective value
* final time
* state trajectories
* control trajectory
* boundary conditions
* path constraints
* terminal conditions
* solver convergence and status

A solver reporting “success” is NOT sufficient. Verify the numerical solution against the MATLAB/reference result within defensible numerical tolerances.

Do not proceed to other examples until this benchmark passes.

# Subsequent validation

After Goddard passes, validate sequentially:

1. `yoptimization/pages/examples/bangFreeTf`
2. `yoptimization/pages/examples/bryson`
3. `yoptimization/pages/examples/goddardRocket`
4. `yoptimization/pages/examples/isoConstraint`
5. `yoptimization/pages/examples/nonDiff`
6. remaining cases under `pages/benchmark`

For every case:

**inspect MATLAB → reproduce the formulation using ****`eco_opt`**** → install/import the package → execute → compare numerically → debug → add regression test → proceed**

Never declare an example working solely because it executes without an exception.

# Testing requirements

Add automated tests covering:

* package installation and import
* public API imports
* variable and expression construction
* dynamics and objective definitions
* bounds and constraints
* discretization/transcription
* solver invocation
* result and trajectory retrieval
* Goddard benchmark regression
* every subsequently validated example

Tests must run with:

```bash
pytest
```

Tests must not depend on MATLAB being installed or available at runtime.

Where exact numerical equality is inappropriate, use documented absolute and relative tolerances.

# Debugging rules

When results disagree:

1. Determine the numerical or mathematical cause before changing code.
2. Compare MATLAB and Python formulations term by term.
3. Check scaling, units, bounds, signs, indexing, interpolation, derivatives, discretization, initial guesses, and solver tolerances.
4. Verify that the installed package is the code being tested.
5. Fix the general implementation rather than introducing example-specific patches.
6. Re-run all previously passing regression tests after significant changes.

# Definition of done

The project is complete only when:

* `eco_opt` is a proper installable Python package;
* `pip install .` succeeds in a clean environment;
* `import eco_opt` succeeds;
* the public API is documented and usable;
* the package does not require MATLAB at runtime;
* the Goddard benchmark reproduces the reference solution;
* all listed examples execute successfully through `eco_opt`;
* numerical outputs agree with their MATLAB/reference counterparts within documented tolerances;
* automated regression tests exist and pass;
* previously passing examples remain passing;
* package metadata and dependencies are declared;
* there are no benchmark-specific hacks masquerading as general functionality;
* the deliverable is the package source tree, not a `.zip` archive of scripts.

Maintain a short progress log showing:

`case | status | objective/reference error | trajectory error | package/test status | remaining issue`

Prioritize a functioning, installable `eco_opt` package and numerical verification over commentary.
