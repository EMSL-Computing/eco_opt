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
`Maximum_Iterations_Exceeded`, …).

## IPOPT options

Pass a dict to `solve` / `optimize` / `closed_loop` as
`ipopt_options`. Keys are IPOPT option names; values are the option
values. CasADi forwards them to IPOPT. The full list is
[Ipopt Options](https://coin-or.github.io/Ipopt/OPTIONS.html)
(Wächter and Biegler, 2006).

```python
sol = ocp.solve(
    control_intervals=40,
    ipopt_options={
        "max_iter": 5000,
        "print_level": 5,
        "tol": 1e-6,
    },
)
```

`eco_opt` sets `print_level` to 0 unless the dict overrides it, and
sets CasADi `print_time` to 0. IPOPT is therefore silent until
`print_level` is raised.

| Option | Role |
|---|---|
| `print_level` | 0 silent … 5 iteration log … 12 verbose |
| `max_iter` | iteration cap (default in IPOPT is 3000) |
| `max_cpu_time` | wall-time cap in seconds |
| `tol` | relative NLP error for `Solve_Succeeded` |
| `acceptable_tol` | looser tolerance for `Solved_To_Acceptable_Level` |
| `constr_viol_tol` | constraint violation |
| `nlp_scaling_method` | `"gradient-based"` (default), `"none"`, `"equilibration-based"` |
| `mu_strategy` | `"monotone"` or `"adaptive"` (barrier) |
| `warm_start_init_point` | `"yes"` / `"no"` — IPOPT dual warm start, not `optimize(warm_start=...)` |

`optimize(warm_start=True)` reuses the **primal** \(x_0\) from the
last solve. That is independent of IPOPT’s
`warm_start_init_point`.

`closed_loop` also defaults `max_iter` to 2000 unless
`ipopt_options` sets it. Mesh defaults for closed loop differ from
`solve`; see [API](api.md).

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
microbial-growth and irrigation examples do that.

## IPOPT missing

If CasADi imports but the solver cannot find IPOPT, the wheel did
not ship a usable `libipopt`. Install IPOPT from the platform
package manager and confirm `import casadi` still works.
