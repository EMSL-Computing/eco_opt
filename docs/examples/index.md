# Examples

Each page uses the same skeleton:

1. Question
2. Formulation (states, control, dynamics, objective, constraints, horizon)
3. What was done
4. What we found
5. Interpretation
6. Figures
7. How to run / tests
8. References — **last**, and only when a paper or Yop page must be cited

Typewriter names match folders under `examples/` in the **git clone**
(not the wheel). Plot scripts need matplotlib; Goddard comparison
plots and Lindh’s 1-D check need scipy:

```bash
python -m pip install -e ".[examples]"
cd examples/microbial_growth
python microbial_growth.py
```

Every **Run and tests** block starts with `cd examples/<folder>` so it
works from the repository root. Plot scripts need `.[examples]`.
Those blocks also run pytest, which is **`[test]`** (pytest plus the
plot extras). Scripts import the public package
(`import eco_opt`, and `from eco_opt import Namespace` when named
signals are used). Do not import `eco_opt.system` or other internal
modules. `Namespace` is names only. Import `casadi as cs` when \(f\)
uses CasADi functions (`exp`, `sin`, `sign`, interpolants); a list of
derivatives needs no CasADi import
([CasADi expressions](../casadi-expressions.md)). Each page has a short public-API
skeleton; the full script is under `examples/<name>/`.

A suggested path: [First solve](../first-solve.md) (one worked script),
[Posing a problem](../posing.md),
[Transcription](../transcription.md),
[Constraints and time](../constraints-and-time.md),
[Build and optimize](../build-optimize.md),
then the tables below. `eco_opt` is not only plants and microbes:
lakes, greenhouse climate, irrigation, harvest, and textbook OCPs
(rockets, bang–bang, path constraints) sit alongside allocation and
uptake models. **Traditional OCP** pages supply a closed form or a
transcription check. [Troubleshooting](../troubleshooting.md) if the
NLP fails.

Yop-shipped problems are rewritten here and **link** the Yop page in
References. The wording is not copied from
[yoptimization.com](https://www.yoptimization.com).

`microbial_growth` is Manzoni, Chakrawal and Ledder (2023), not a Yop
website example. Cases inspired by other papers state in References that
trajectories are not validated against those papers’ figures
(Chakrawal et al. excepted).

## Ecological models

| Example | Question |
|---|---|
| [microbial_growth](microbial_growth.md) | When should uptake stop as substrate is exhausted? |
| [opt_lignin](opt_lignin.md) | When to pay for ligninolytic oxidation? |
| [root_opt](root_opt.md) | How to split belowground carbon under N limitation? |
| [cohen_king_roughgarden](cohen_king_roughgarden.md) | Vegetative vs reproductive allocation |
| [iwasa_roughgarden_shoot_root](iwasa_roughgarden_shoot_root.md) | Shoot, root, and fruit |
| [lindh_2016_flowering](lindh_2016_flowering.md) | Flowering time under constrained growth |
| [multiphase_cohen](multiphase_cohen.md) | Cohen’s switch as two phases |
| [shallow_lake](shallow_lake.md) | Phosphorus loading and two management locals |
| [greenhouse](greenhouse.md) | Heating with tabulated weather |
| [irrigation_bucket](irrigation_bucket.md) | Irrigation when rain is not in the forecast |

## Traditional OCP problems

Textbook and benchmark transcriptions. They are not ecological models.
They check free time, path constraints, bang–bang structure, and
collocation against closed forms or an archive.

| Example | Question |
|---|---|
| [goddard_benchmark](goddard_benchmark.md) | Dimensional Goddard, free horizon |
| [goddard_rocket_basic](goddard_rocket_basic.md) | Nondimensional Goddard |
| [multiphase_goddard_landing](multiphase_goddard_landing.md) | Ascent then landing |
| [bang_free_tf](bang_free_tf.md) | Minimum-time double integrator |
| [bryson_denham](bryson_denham.md) | Path-constrained minimum energy |
| [iso_constraint](iso_constraint.md) | Fixed integral of \(u^2\) |
| [non_diff](non_diff.md) | Discontinuous pulse in \(f\) |
| [brachistochrone](brachistochrone.md) | Minimum-time cycloid |
| [van_der_pol](van_der_pol.md) | Forced oscillator, running cost |
| [lenhart_workman](lenhart_workman.md) | Three textbook OCPs |
| [lenhart_workman_labs](lenhart_workman_labs.md) | Book labs, posed in eco_opt |
| [xuezhong_wang](xuezhong_wang.md) | Fixed-time double integrator |

## Numerical features

| Example | Question |
|---|---|
| [algebraic_qss](algebraic_qss.md) | Fast biomass as an algebraic closure |

Multi-phase collocation is shown inside
[multiphase_cohen](multiphase_cohen.md) and
[Goddard landing](multiphase_goddard_landing.md). Open vs closed loop
is shown inside the irrigation, microbial, lignin, and root-opt
pages.
