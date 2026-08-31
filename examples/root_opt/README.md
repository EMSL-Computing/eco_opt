# root_opt

Python/`eco_opt` port of the optimality-based plant carbon-nitrogen
allocation model from Chakrawal, Mooney & Colombi (2025/2026), "Modelling
belowground plant acclimation to low soil nitrogen — an eco-evolutionary
approach", *New Phytologist* 251:555-568,
https://doi.org/10.1111/nph.71210 (code: `external/root_opt/`,
https://doi.org/10.5281/zenodo.17872549).

## Problem

    maximize   int_0^T leaf_growth_rate dt
    subject to  9-state plant-soil-microbe C/N ODE system (leafC, rootC,
                exuC, SOMC, micC, rootN, SOMN, micN, inorgN)
                2 controls: rootCSupply [gC/m2/d], exu_fraction [-]
                (root exudation modeled as a *fraction* of net root C
                assimilation -- a nested/hierarchical control structure
                that couples exudation and root growth and avoids
                artificial bang-bang competition between the two sinks)
                various state/path bounds (C:N ratio caps, minimum
                leaf/root growth allocation fractions, N-uptake-capacity
                normalization, exudation-fraction bounds [0.1, 0.5])
                x(0) = base initial condition (leafC=50, rootC=75,
                exuC=0.5, SOMC=7200, micC=2.8, inorgN=20 gN/m2, plus
                stoichiometric N pools), fixed t_f = 120 d

Port of `external/root_opt/solve_ocp_nested.m` +
`params_base_condition.m` (the "nested"/Option-1 control structure from
the manuscript). `plant_soil_nitrogen_allocation.py` implements the
dynamics 1:1 (same equations, parameter values, and bounds) and
`build_ocp`/`solve` construct/solve the corresponding eco_opt `Ocp`.

## Validation

No closed-form reference exists for this 9-state ecological OCP (as in
the MATLAB original). We validate:

1. **Convergence**: `Solve_Succeeded` for the base scenario and two
   N-stress scenarios (low mineral N, high SOM C:N).
2. **Regression**: base-scenario objective tracked at
   `REGRESSION_OBJECTIVE = 499.505` (control_intervals=50, degree=3,
   legendre collocation).
3. **Mass-balance / physical sanity**: all 9 state pools stay
   non-negative; leaf C and root C grow substantially over the season;
   mineral N is net depleted (uptake exceeds mineralization, consistent
   with `leaching_rate=0`); control bounds respected.
4. **Ecological consistency (N-limitation)**: reproducing
   `Figure2_and_3.m`'s three-scenario comparison (base / low mineral N /
   high SOM C:N) in `plot_root_opt.py` shows the physically-expected
   pattern -- N-limited scenarios have *lower* cumulative leaf growth and
   *higher* exudation-fraction control (investing more in N acquisition
   under stress) than the base scenario. See `figs/root_opt_scenarios.png`.
   `tests/test_root_opt.py::test_low_mineral_n_scenario_reduces_growth`
   pins this qualitative direction as a regression guard.

## Running

```bash
python plant_soil_nitrogen_allocation.py   # base scenario, prints summary
python plot_root_opt.py                    # writes figs/root_opt_scenarios.png (custom summary figure)
python reproduce_matlab_figures.py         # writes figs/Figure2_revised.png and figs/Figure3.png,
                                            # matching the *exact* layout/labels/colors of
                                            # external/root_opt/Figure2_and_3.m's exported figures
pytest ../../tests/test_root_opt.py -q
```

`reproduce_matlab_figures.py` mirrors `Figure2_and_3.m` panel-for-panel:
same 5x3 (Figure2_revised) / 4x1 (Figure3) tiled layout, panel order and
(A)-(J)/(A)-(D) labels, `copper(3)` 3-scenario color palette, axis
labels (Anet, G_L/A_net, S/A_net, G_R/S_net, E/S_net time series +
cumulative bars; root N uptake/growth/exudation rate/net mineralisation
in Figure3), dashed reference lines (0.2/0.8 thresholds,
min/max_exudation bounds), hidden x-tick labels on non-bottom rows, and
the same solver settings as the MATLAB script (`control_intervals=100`).
`plot_root_opt.py` is a separate, simpler custom summary figure (states +
controls + objective bar) kept for quick sanity-checking during
development.


## Notes / simplifications vs. the MATLAB original

* This port implements the base/N-stress *scenario comparison* (as in
  `Figure2_and_3.m`'s first figure) but not the full parameter-sweep
  figures (`Figure_4/5/6_*.m`, which vary soil C:N x initial N x plant
  traits over dozens of combinations using pre-computed `.mat` result
  caches) -- those would be a natural follow-up but are out of scope for
  a single example script.
* No warm-start / mesh continuation is implemented (the MATLAB script's
  `tdf` warm-start branch in `solve_ocp_nested.m`); each scenario here is
  solved cold-start at `control_intervals=50, polynomial_degree=3`, which
  converges directly for all three scenarios shown.
