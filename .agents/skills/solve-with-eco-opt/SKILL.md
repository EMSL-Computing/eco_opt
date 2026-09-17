---
name: solve-with-eco-opt
description: Formulate, solve, simulate, troubleshoot, and interpret ecological optimal-control models using eco_opt. Use when a researcher provides equations, states, controls, parameters, constraints, data, or a scientific optimization question. Do not use for maintaining eco_opt internals.
---

# Solve with eco_opt

Help the user turn an ecological or process-based model into a correct,
runnable eco_opt workflow.

## Establish the problem

Determine from the information available:

- states and their initial conditions;
- controls and their bounds;
- model parameters and units;
- state dynamics;
- algebraic variables or equations, if any;
- path and endpoint constraints;
- objective and optimization direction;
- fixed or free final time; and
- desired outputs, simulations, or plots.

Ask for missing scientific information only when it materially changes the
mathematical problem. Otherwise, state reasonable assumptions.

## Build the model

Use only the public API imported through `import eco_opt`. Follow the public API
and repository-wide requirements in `AGENTS.md`.

Return model derivatives as a scalar, list, or tuple. Do not use
`casadi.vertcat` for derivatives in user-facing code. Use CasADi where symbolic
functions such as `exp`, `sin`, `sign`, or interpolants are needed.

Prefer a complete, executable example over isolated fragments. Use the
smallest model structure that expresses the scientific problem, and explain
the mapping from the equations to the code.

## Read relevant repository guidance

Read only the documentation needed for the current task, resolving these paths
from the repository root:

- Basic formulation: `docs/posing.md`
- Symbolic expressions: `docs/casadi-expressions.md`
- Constraints and time: `docs/constraints-and-time.md`
- First solution: `docs/first-solve.md`
- Build, constraint, and optimization sweeps: `docs/build-optimize.md`
- Initial guesses, mesh, and scaling: `docs/guesses-mesh-scaling.md`
- Results and plots: `docs/results-and-plots.md`
- Solver failures: `docs/troubleshooting.md`
- API details: `docs/api.md`

Inspect an existing example with similar mathematical features when that is
more useful than starting from an empty template.

When plots would help the user understand or assess a solution, read and apply
the [standard plotting guide](references/plotting.md). Adapt labels and units
to the scientific model rather than exposing anonymous array indices.

## Present the result

Explain:

1. how the scientific model maps to states, controls, and constraints;
2. any assumptions introduced;
3. how to run the code;
4. what quantities should be inspected to assess the solution; and
5. any numerical limitations or convergence concerns.

Do not claim scientific validity merely because the optimizer converged.
Distinguish numerical success from biological interpretation.

Do not modify the eco_opt library unless the user explicitly requests a
library change.
