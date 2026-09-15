# eco_opt

`eco_opt` is a Python package for turning a process-based ODE into a
constrained optimal-control problem and solving it by direct collocation.

The process right-hand side is supplied by the user. Controls
(quantities the modeller lets vary in time to improve \(J\)), the
performance criterion, and constraints are declared. The package
transcribes that statement with CasADi and solves it with IPOPT.

The modelling workflow follows [Yop](https://www.yoptimization.com)
(Leek, 2016). `eco_opt` is a Python implementation of that workflow, not
a copy of the Yop documentation. When an example here is a transcription
of a Yop-shipped case, the example page links to the corresponding Yop
page.

## Start here

1. [Install](install.md)
2. [First solve](first-solve.md) — microbial growth (finite substrate)
3. [Results and plots](results-and-plots.md) — `sol` arrays and matplotlib
4. [Posing a problem](posing.md) — states, controls, costs
5. [CasADi expressions](casadi-expressions.md) — `cs.exp`, `if_else`, interpolants in the model
6. [Transcription](transcription.md) — piecewise-constant \(u\), collocation, mesh
7. [Constraints and time](constraints-and-time.md) — boxes, path, \(t_0\)/\(t_f\), horizon
8. [Guesses, mesh, scaling](guesses-mesh-scaling.md) — cold start, mesh, units
9. [Build and optimize](build-optimize.md) — frozen mesh, bound and parameter sweeps
10. [Examples](examples/index.md) — ecological models, then textbook OCPs
11. [Troubleshooting](troubleshooting.md)
12. [API](api.md)

## What it is for

Process models in ecology are usually integrated forward: given the
process equations and inputs, stocks (state variables) are predicted.
Some of those inputs or rates can be treated as controls. `eco_opt`
asks what trajectory is implied by an explicit optimality hypothesis
and the same process equations, subject to constraints.

The same workflow applies to populations, water, greenhouse climate,
biogeochemistry, and other process-based ecological ODEs — not only
carbon allocation or organism physiology.

It does not claim that the system maximises the chosen objective. A
numerical solution shows that the hypothesis is internally consistent
with the dynamics and constraints.

## Cite

See [Cite](cite.md). Also cite CasADi, IPOPT, and Yop for the numerical
stack and modelling workflow.
