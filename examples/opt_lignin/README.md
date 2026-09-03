# Ligninolytic investment

Optimal ligninolytic activity during litter decomposition (Chakrawal et
al. 2024, https://doi.org/10.1111/nph.19572).

## Problem

Maximize cumulative microbial growth `G` over a fixed horizon `TT`.
Litter carbon is split into total pool `CT` and aromatic sub-pool `Co`;
the control `vo` is investment in lignin-targeting decomposition.

```
maximize   int_0^TT G dt
subject to CTdot = -Dh - Do + G,   CT(0) = CT_0
           Codot = -Do + mo*G,     Co(0) = Co_0 = CT_0 * ARC0
           0 <= CT, 0 <= Co
           0 <= vo <= vomax
           t_f == TT

L  = Co / CT
Dh = vh_max * g(L) * (CT - Co)
Do = vo * Co
e  = emax - ro * vo
G  = e * (Dh + Do)
g(L) = exp(-(L/a)^b)
```

`vh_max`, `mo`, `ro` are fixed parameters (pinned via box-equality
constraints on `sys.p`). `plot_figure2.py` compiles the NLP once and
retunes those pins with `set_constraints` / `optimize`.

## Run

```bash
python lignin_decomposition.py
python plot_figure2.py
```

## Validation

No closed-form objective is published. Checks: solver success across the
parameter sweep; non-negative pools; delayed onset of `vo`;
regression-tracked baseline objective.

## Baseline parameters

| vh_max | mo  | ro | ARC0 | CT_0 | TT   | CN0 |
|--------|-----|----|------|------|------|-----|
| 0.0015 | 0.1 | 25 | 0.05 | 0.5  | 2000 | 50  |
