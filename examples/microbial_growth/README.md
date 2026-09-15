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
cd examples/microbial_growth
python microbial_growth.py
python plot_microbial_growth.py
python open_vs_closed.py
python -m pytest ../../tests/test_microbial_growth.py -q
```

## Validation

Collocation is checked against an independent PMP trajectory
(`analytical_solution`; figure `figs/microbial_growth_vs_analytical.png`)
and against resource exhaustion / box bounds. $T^*\approx 1.518$ and
$J^*\approx 0.1518$ are also regression-tracked.

## Open loop vs closed loop

`open_vs_closed.py` pins $T=1.5$ (the free-$T$ optimum is $\approx 1.52$)
and compares:

1. One open-loop $u^*(t)$, replayed on the same $f$.
2. Shrinking-horizon closed loop on the same $f$.
3. The same $u^*(t)$ applied to a plant that gets extra C on
   $[0.4,0.7]$, versus re-planning from the measured $C$.

Same plant: the two $u$ series match and $\int G=0.152$. Pulsed plant:
the tape still uses the planned uptake and leaves $C(T)\approx 0.22$
uneaten ($\int G=0.152$); closed loop raises $u$ after the pulse,
exhausts C, and reaches $\int G=0.190$.

Figure: `figs/open_vs_closed.png`.

## Parameters (baseline `k=0`)

| beta | gamma | r   | emax | mu  | k |
|------|-------|-----|------|-----|---|
| 0.5  | 0.5   | 0.1 | 0.5  | 0.5 | 0 |
