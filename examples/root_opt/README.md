# Plant–soil C–N allocation

Optimality-based plant carbon–nitrogen allocation (Chakrawal, Mooney &
Colombi 2026, *New Phytologist* 249:555–568,
https://doi.org/10.1111/nph.71210).

## Problem

```
maximize   int_0^T leaf_growth_rate dt
subject to  9-state plant-soil-microbe C/N ODE system
            2 controls: rootCSupply [gC/m2/d], exu_fraction [-]
            exudation as a fraction of net root C assimilation
            t_f = 120 d
```

## Validation

No closed-form reference exists. Checks: solver success for base and two
N-stress scenarios; regression-tracked base objective; non-negative
pools; N-limited scenarios have lower cumulative leaf growth and higher
exudation fraction.

## Run

```bash
python plant_soil_nitrogen_allocation.py
python plot_root_opt.py
python reproduce_figures.py
pytest ../../tests/test_root_opt.py -q
```
