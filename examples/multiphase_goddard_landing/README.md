# Multi-phase Goddard landing

Two-phase Goddard rocket: climb to a free apogee, then land softly.
This is the Yop example
[Goddard Rocket Landing](https://www.yoptimization.com/goddardLanding)
transcribed with `eco_opt.Phase`. Same dimensional rocket as
`examples/goddard_benchmark/` (max ascent, no landing).

---

## 1. Problem

The rocket is the Bryson–Ho / Yop Goddard vehicle. States are velocity
$v$, altitude $h$, mass $m$; control is fuel mass-flow $F$:

$$
\dot v = \frac{F c - \mathrm{sign}(v)\,D_0 e^{-\beta h} v^2}{m} - g(h),
\quad
\dot h = v,
\quad
\dot m = -F,
\quad
g(h)=g_0\bigl(r_0/(r_0+h)\bigr)^2.
$$

Unlike the single-phase max-ascent benchmark, the vehicle must return
to the pad:

| Phase | Objective | Terminal |
|---|---|---|
| 1 ascent | maximize $h(t_s)$ | $t_s$ free |
| 2 landing | none (Yop writes $\min t_f(0)$) | $h(T)=0$, $v(T)=0$, $T$ free |

Junction: $t_s$ shared, $(v,h,m)$ continuous. Path bounds on both
phases: $h\ge 0$, $m_f \le m \le m_0$, $0\le F\le F_m$.

The markdown problem statement also lists $m(T)=m_f$. The Yop
**code** does not pin terminal mass, only the path interval. This
example follows the code: leftover fuel is allowed. At the default
mesh the optimizer still burns down to $m_f$.

$\mathrm{sign}(v)$ in the drag term is required on the descent
($v<0$); the max-ascent benchmark never goes negative.

---

## 2. Why two phases

A single grid with $h(T)=v(T)=0$ and $\max h(t)$ is an
*interior* peak, which collocation can represent, but the switch
between “go up” and “come down” is a change of *task*, not just of
$u$. Splitting at a free $t_s$ and maximising $h(t_s)$ makes the
apogee a decision variable and gives each arc its own collocation
grid. That is the Yop listing (`ocp(1)`, `ocp(2)`).

---

## 3. `eco_opt` transcription

Two `System` instances (distinct symbols) share the same model
callback. The Mayer term is attached to **phase 1** because it is
built from `sys1.y.rocket.height`:

```python
sys1 = System(..., name="goddard_ascent")
sys2 = System(..., name="goddard_landing")  # same model, new symbols
p1 = Phase(sys1, constraints=[t_0 conditions, path bounds])
p2 = Phase(sys2, constraints=[path bounds, t_f(v)==0, t_f(h)==0])
ocp.max(t_f(sys1.y.rocket.height))
ocp.st(phases=[p1, p2], continuity=True)
sol = ocp.solve(collocation_points="radau", tf_guess=[0, 40, 80])
# or: ocp.build(..., tf_guess=[0, 40, 80]); sol = ocp.optimize()
```

Yop scales the objective by $10^{-4}$. `eco_opt` has no scaling API
yet; IPOPT `nlp_scaling_method=gradient-based` is enough on this
mesh.

---

## 4. Default mesh and regression

Radau, degree 3, 50 intervals per phase (Yop: 100 intervals, same
degree, Radau). Verified solve:

| quantity | value |
|---|---|
| peak $h(t_s)$ | $8.208\times 10^{4}$ m |
| $t_s$ | 153.5 s |
| $T$ | 289.3 s |
| $h(T), v(T)$ | 0 |
| $m(T)$ | $m_f$ |

For comparison, the single-phase max-ascent benchmark (no landing,
$v\ge 0$) reaches $h(T)\approx 1.61\times 10^{5}$ m. The landing
problem is lower: fuel is reserved for the descent.

No published closed form. Tests regression-track peak height and $T$
at this mesh, plus the landing equalities and bounds.

---

## 5. Files

| File | Role |
|---|---|
| `goddard_landing.py` | two-phase OCP |
| `plot_goddard_landing.py` | states/control and drag–$h$–$v$ 3-D plot |
| `tests/test_multiphase_goddard_landing.py` | landing + regression |

---

## 6. Running

```bash
cd examples/multiphase_goddard_landing
python goddard_landing.py
python plot_goddard_landing.py
python -m pytest ../../tests/test_multiphase_goddard_landing.py -q
```

---

## 7. Parameters

| symbol | value | meaning |
|---|---|---|
| $m_0$ | 214.839 kg | initial mass |
| $m_f$ | 67.9833 kg | dry mass (lower bound) |
| $F_m$ | 9.525515 kg/s | max fuel flow |
| $c$ | 2060 m/s | exhaust speed |
| $D_0$ | 0.01227 | drag coefficient |
| $\beta$ | $0.145\times 10^{-3}$ 1/m | density scale |
| $g_0$ | 9.81 m/s$^2$ | sea-level gravity |
| $r_0$ | $6.371\times 10^6$ m | Earth radius |

---

## 8. References

- Bryson, A. E. & Ho, Y.-C. (1975). *Applied Optimal Control*.
- Yop example: Goddard rocket landing
  (https://www.yoptimization.com/goddardLanding).
- Single-phase max-ascent sibling: `examples/goddard_benchmark/`.
