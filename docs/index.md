# eco_opt

`eco_opt` is a Python package for turning a process-based ODE into a
constrained optimal-control problem and solving it by direct collocation.

You keep the mass-balance right-hand side you already write. You declare
which fluxes may vary in time, which performance criterion is being
tested, and which physiological limits close the problem. The package
transcribes that statement with CasADi and solves it with IPOPT.

The modelling workflow follows [Yop](https://www.yoptimization.com)
(Leek, 2016). `eco_opt` is a Python implementation of that workflow, not
a copy of the Yop documentation. When an example here is a transcription
of a Yop-shipped case, the example page links to the corresponding Yop
page.

## Start here

1. [Install](install.md)
2. [First solve](first-solve.md) — Goddard rocket (traditional OCP)
3. [Posing a problem](posing.md) — states, controls, costs
4. [Examples](examples/index.md) — ecological models, then textbook OCPs
5. [Troubleshooting](troubleshooting.md)
6. [API](api.md)

## What it is for

Process models in ecology are usually integrated forward: rates and
allocation fractions are specified, pools (states) are predicted. Many of
those fluxes can vary biologically. `eco_opt` asks what trajectory is
implied by an explicit optimality hypothesis and the same process
equations.

It does not claim that organisms maximise the chosen objective. A
numerical solution shows that the hypothesis is internally consistent
with the dynamics and constraints.

## Cite

See [Cite](cite.md). Also cite CasADi, IPOPT, and Yop for the numerical
stack and modelling workflow.
