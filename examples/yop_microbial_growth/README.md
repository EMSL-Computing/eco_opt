# Microbial growth

Finite-substrate microbial growth (Manzoni, Chakrawal & Ledder 2023,
https://doi.org/10.3389/fevo.2023.1094269).

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

`C` is the substrate pool (state), `u` the uptake rate (control). Free
terminal time: the optimizer finds when continued uptake is no longer
worthwhile as the pool nears exhaustion.

## Run

```bash
/home/chak803/miniforge3/envs/llm_env/bin/python microbial_growth.py
/home/chak803/miniforge3/envs/llm_env/bin/python plot_microbial_growth.py
```

## Validation

No closed-form objective is published. Checks: solver success; `C(T) -> 0`;
bounds respected; regression-tracked objective and terminal time.

## Parameters (baseline `k=0`)

| beta | gamma | r   | emax | mu  | k |
|------|-------|-----|------|-----|---|
| 0.5  | 0.5   | 0.1 | 0.5  | 0.5 | 0 |
