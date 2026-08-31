# opt_lignin_figure2

Python/`eco_opt` port of the OptLignin "Figure 2" optimal-control model from
Chakrawal et al. 2024, "Modelling optimal ligninolytic activity during plant
litter decomposition" (https://doi.org/10.5281/zenodo.10553503), originally:

```
examples/OptLignin/OptLignin1.0/
    opt_con.m       -- YopSystem/YopOcp model (fixed final time TT)
    Figure2.m        -- 2x4 parameter-sweep figure (vh_max, ro, ARC0, mo)
    set_up.m         -- shared parameters / g-function definition
```

## Problem

Maximize cumulative microbial growth `G` over a fixed decomposition horizon
`TT`, where litter carbon is split into a total pool `CT` and an aromatic
("lignin-like") sub-pool `Co`, and the control `vo` is the rate of
investment in lignin-targeting decomposition:

```
maximize   int_0^TT G dt
subject to CTdot = -Dh - Do + G,   CT(0) = CT_0
           Codot = -Do + mo*G,     Co(0) = Co_0 = CT_0 * ARC0
           0 <= CT, 0 <= Co
           0 <= vo <= vomax
           t_f == TT  (fixed final time)

L  = Co / CT
Dh = vh_max * g(L) * (CT - Co)      # non-aromatic decomposition
Do = vo * Co                        # aromatic-targeting decomposition (control-driven)
e  = emax - ro * vo                 # growth efficiency (decreasing in effort)
G  = e * (Dh + Do)
g(L) = exp(-(L/a)^b)                # decomposability-vs-aromaticity kernel
```

`vh_max`, `mo`, `ro` are treated as fixed parameters (Yop's
`ocp.Parameter(...).setBound(...)` pattern -- pinned via box-equality
constraints on `sys.p`, exactly mirroring `opt_con.m`).

## Package fix required

Porting this example surfaced a real gap in `eco_opt.System`: models that
declare `parameters=N` but no algebraic states (`algebraics=0`) never
received the parameter vector `p` in their model callback (only the
`algebraics>0` branch passed `p`). Fixed in `src/eco_opt/system.py` so any
model with `parameters>0` gets called as `model(t, x, u, p)`.

## Files

* `lignin_decomposition.py` -- `eco_opt` port (`System`/`Ocp`/`solve()`), 30
  control intervals, degree-3 legendre collocation (matches `opt_con.m`'s
  `'controlIntervals', param.numiter (=30), 'polynomialDegree', 3,
  'collocationPoints', 'legendre'`).
* `plot_figure2.py` -- reproduces a Python analogue of the paper's Figure 2:
  a 2x4 grid sweeping vh_max, ro, ARC0, mo (each over two values), plotting
  CT/CT(0), Co/Co(0) (top row) and vo(t) (bottom row).
* `../../tests/test_opt_lignin_figure2.py` -- regression test.

## Run

```bash
/home/chak803/miniforge3/envs/llm_env/bin/python lignin_decomposition.py
/home/chak803/miniforge3/envs/llm_env/bin/python plot_figure2.py
```

Output figure: `figs/opt_lignin_figure2.png`.

## Validation

No single closed-form objective is published; the paper's Figure 2 shows
qualitative sensitivity trends only. We validate:

* **success**: IPOPT reports `Solve_Succeeded` for the baseline configuration
  and for all 8 parameter-sweep solves in `plot_figure2.py`.
* **mass-balance sanity**: `CT >= 0`, `Co >= 0` throughout every trajectory
  (checked in the test).
* **qualitative match**: control `vo` stays near zero early in the horizon
  and ramps up mid-decomposition before tapering near `t=TT`, consistent
  with the paper's framing of delayed onset of ligninolytic activity
  (spot-checked visually against `figs/opt_lignin_figure2.png`).
* **regression tracking**: `REGRESSION_OBJECTIVE = 0.16281350935527988` at
  the baseline parameters (`vh_max=0.0015, mo=0.1, ro=25, CT_0=0.5,
  ARC0=0.05, TT=2000`) guards against silent future regressions.

## Baseline parameters (Figure2.m defaults)

| vh_max | mo  | ro | ARC0 | CT_0 | TT   | CN0 |
|--------|-----|----|------|------|------|-----|
| 0.0015 | 0.1 | 25 | 0.05 | 0.5  | 2000 | 50  |
