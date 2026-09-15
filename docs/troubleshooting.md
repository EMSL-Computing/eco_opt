# Troubleshooting

Direct collocation is a nonlinear program. `sol.success` false, a
restoration failure, or a trajectory that ignores a bound usually
means the **NLP** is poorly scaled or poorly guessed, not that the
hypothesis is empty.

Mesh, `tf_guess`, constant `initial_guess`, warm starts, and
nondimensional units:
[Guesses, mesh, scaling](guesses-mesh-scaling.md).
Piecewise-constant \(u\) and CasADi `if_else`:
[Transcription](transcription.md).
Boxes vs path vs `t_0`/`t_f`, and explicit time:
[Constraints and time](constraints-and-time.md).

## Read the solver first

```python
print(sol.success, sol.return_status, sol.tf, sol.objective)
```

`return_status` is IPOPT’s string (`Solve_Succeeded`,
`Solved_To_Acceptable_Level`, `Restoration_Failed`,
`Maximum_Iterations_Exceeded`, …). Pass extra keys in
`ipopt_options`, for example `{"max_iter": 5000, "print_level": 5}`.

## Restoration failed / iteration limit

Useful sequence:

1. Coarser mesh, then refine (new `solve`; same constant guess, a
   `tf_guess` near the coarse \(T^*\)).
2. A `tf_guess` of the correct order of magnitude and box bounds that
   contain the guess.
3. Radau when \(u^*\) is expected to switch.
4. `print_level=5` to identify an active constraint.
5. Consistency of an algebraic residual \(g=0\) (`algebraics`; an
   initial \(z(0)\) only when it is part of the formulation).

On a **frozen** mesh, `build` / `set_constraints` / `optimize`
warm-starts a bound or parameter change. A new mesh rebuilds the NLP.

## Closed loop

If `closed_loop` disagrees with `solve` on the **same** process, the
sampling interval `dt` is too large or the remaining horizon no
longer contains the planned terminal time (free \(T\) with a
shrinking window past burnout). A fixed \(T\) is required for the
closed-loop study if the open-loop problem had a free horizon; the
microbial example does that.

## IPOPT missing

If CasADi imports but the solver cannot find IPOPT, the wheel did
not ship a usable `libipopt`. Install IPOPT from the platform
package manager and confirm `import casadi` still works.
