# yop_microbial_growth

Python/`eco_opt` port of the MATLAB Yop microbial-growth optimal-control
example from Manzoni et al. 2022 (Frontiers in Ecology and Evolution,
https://doi.org/10.3389/fevo.2023.1094269), originally located at:

```
examples/OptLignin/YOP benchmark/Manzoni et al. 2022/
    GrowthModel.m, yop_opt.m, YOP_microbial_growth.m   (Yop/CasADi direct)
    indirect_bvp4c.m                                    (indirect shooting, bvp4c)
    indirect_forward_backward_sweep.m                   (Lenhart & Workman FBS)
```

## Problem

Maximize cumulative microbial growth `G` over a free time horizon `[0, T]`,
subject to first-order resource-pool dynamics with decay and necromass
recycling:

```
maximize   int_0^T G(u) dt
subject to Cdot = -u - gamma*C + mu*G(u),   C(0) = 1,   0 <= C <= 1
           G(u) = emax*beta*(u - r)/(u + beta) - k*u
           0 <= u <= 10
           T free, T >= 0
```

`C` is the substrate pool (state), `u` the uptake rate (control), `G` the
growth-rate return on uptake. Free terminal time: the optimizer implicitly
finds the time at which continued uptake is no longer worthwhile as the
resource pool nears exhaustion.

## Files

* `microbial_growth.py` -- `eco_opt` port (`System`/`Ocp`/`solve()`), 200
  control intervals, radau collocation (matches `yop_opt.m`'s
  `'controlIntervals', 200, 'collocationPoints', 'radau'`).
* `plot_microbial_growth.py` -- generates comparison figures analogous to
  the MATLAB `indirect_bvp4c_Fig1.png` / `indirect_bvp4c_Fig2.png` plots,
  from eco_opt's own converged solution.
* `../../tests/test_yop_microbial_growth.py` -- regression test.

## Run

```bash
/home/chak803/miniforge3/envs/llm_env/bin/python microbial_growth.py
/home/chak803/miniforge3/envs/llm_env/bin/python plot_microbial_growth.py
```

Output figures land in `figs/microbial_growth_states_control.png` and
`figs/microbial_growth_phase_portrait.png`.

## Validation

No closed-form objective value is published for this exact problem
instance (the paper cross-validates three independent solvers against each
other rather than an analytic reference). We validate:

* **success**: IPOPT reports `Solve_Succeeded`.
* **physical sense**: the free terminal time is driven by resource
  exhaustion, `C(T) -> 0` (checked to `< 1e-4`) -- consistent with the
  paper's framing of the problem as "microbial growth optimization for
  finite food source".
* **box feasibility**: `0 <= u <= 10` and `0 <= C <= 1` are respected
  throughout the trajectory.
* **regression tracking**: `REGRESSION_OBJECTIVE = 0.15184969799373915`,
  `REGRESSION_TF = 1.5181904140431468` (computed with eco_opt at the above
  configuration) guard against silent future regressions in the
  transcription.

## Parameters (baseline `k=0` case, from `YOP_microbial_growth.m`)

| beta | gamma | r   | emax | mu  | k |
|------|-------|-----|------|-----|---|
| 0.5  | 0.5   | 0.1 | 0.5  | 0.5 | 0 |
