# lindh_2016_flowering

Python/`eco_opt` reproduction of Lindh, Johansson, Bolmgren, Lundstrom,
Brannstrom & Jonzen (2016), *Constrained growth flips the direction of
optimal phenological responses among annual plants*, New Phytologist
209:1591-1599, doi:10.1111/nph.13706.

## Problem

Optimal allocation of biomass production `F(V)` between vegetative
growth `V(t)` and reproductive growth `R(t)` over a fixed season
`[0, T]`, via control `u(t) in [0,1]`:

    dV/dt = u(t) F(V),       V(0) = V0
    dR/dt = (1-u(t)) F(V),   R(0) = 0
    maximize R(T)

for three growth models `F(V) = P * U(V)`:

  * exponential (unconstrained): `U(V) = V`
  * logistic (constrained): `U(V) = V(1 - V/Vmax)`
  * WBE (constrained, metabolic scaling): `U(V) = V^0.75 (1 - (V/Vmax)^0.25)`

## Key insight reproduced

The paper proves the optimal control is **bang-bang** (Theorem 1): pure
vegetative growth until a single "flowering time" `tF`, then pure
reproduction. Given that, the central qualitative result (Fig. 2) is
that the *direction* of the optimal flowering-time response to
increasing productivity `P` depends on whether growth is constrained:

  * exponential growth: `tF*(P)` is monotonically **increasing**
  * logistic/WBE (constrained) growth: `tF*(P)` is **hump-shaped**
    (increases, then decreases past a threshold productivity)

i.e. constrained-growth plants eventually flower *earlier* as
productivity keeps rising, while unconstrained-growth plants always
flower later -- the paper's namesake "flip".

## Files

* `optimal_flowering.py`: builds/solves the continuous-relaxation OCP
  (`u(t) in [0,1]`, NOT pre-restricted to bang-bang) via `eco_opt`, for
  all three growth models; verifies numerically that the optimizer's
  `u(t)` is bang-bang (saturates at 0/1, single switch) -- confirming
  Theorem 1 rather than assuming it. Also provides a fast closed-form/
  1-D cross-check: given bang-bang control, `V(t)` obeys the
  *uncontrolled* ODE `dV/dt=F(V)` up to `tF`, so
  `tF* = argmax_tF F(V(tF))(T - tF)` reduces to scalar optimization
  (used to sweep `P` densely in the figure).
* `plot_figure2.py`: reproduces Fig. 2 panels (a)-(i): relative growth
  rate vs `V` (a-c), vegetative-mass trajectories with flowering-time
  markers (d-f), and **optimal flowering time vs productivity** (g-i,
  the paper's central result) for the three growth models -- writes
  `figs/figure2_reproduction.png`.
* `tests/test_lindh_2016_flowering.py`: (1) bang-bang verification for
  all three models via the full OCP; (2) OCP vs 1-D cross-check
  agreement; (3) exponential growth's `tF*(P)` is monotone increasing;
  (4) logistic/WBE growth's `tF*(P)` has an interior maximum and
  decreases afterward (the flip).

## Running

```bash
python optimal_flowering.py    # OCP solves + bang-bang check for P=2, all 3 models
python plot_figure2.py         # writes figs/figure2_reproduction.png
pytest ../../tests/test_lindh_2016_flowering.py -q
```

## Parameters

`Vmax=1`, `V0=0.02`, season `[0,1]` (matches Fig. 2 caption). Productivity
`P` is swept over `[0.3, 15]` in the figure/tests -- wider than the
paper's implicit range, needed because the WBE model's interior maximum
in `tF*(P)` occurs at higher `P` (~4) than logistic's (~3.8) with these
parameters, and a wide sweep is needed to see the post-peak decrease
clearly for both.

## Notes / simplifications

* Only Scenario 1 (Fig. 2: fixed season length, varying productivity)
  is reproduced. Scenario 2 (Fig. 3: varying season length jointly
  with productivity) and the constant-mortality extension (Methods S1 /
  Fig. S1) are out of scope for this example.
* The paper's Methods S2 closed-form expressions for `tF*(P)` are not
  in the extracted PDF text; we instead validate numerically via two
  independent methods (full OCP + 1-D scalar optimization) that agree
  to <0.05 in `tF`.
