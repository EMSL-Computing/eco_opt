# Multi-phase Cohen (1971)

Cohen’s vegetative → reproductive switch transcribed as **two
sequential phases** instead of one collocation grid with a jumping
control. This is the first `eco_opt` example that uses
`eco_opt.Phase` and `Ocp.st(phases=[...])`.

The single-phase sibling (control $u(t)\in[0,1]$ on one horizon) is
`examples/cohen_king_roughgarden/`. Same plant, same closed form; the
point here is the transcription, not the biology.

---

## 1. What “multi-phase” means here

A single-phase OCP has one right-hand side $f$ on $[t_0,t_f]$. A
multi-phase OCP is several of those in sequence, glued by junction
conditions. Cohen’s optimum is bang-bang with a single switch, so it
splits cleanly:

| Phase | Time | Dynamics (bangs baked in) | Control |
|---|---|---|---|
| 1 vegetative | $[0, t_s]$ | $\dot x_1=x_1$, $\dot x_2=0$ | none |
| 2 reproductive | $[t_s, T]$ | $\dot x_1=0$, $\dot x_2=x_1$ | none |

$t_s$ is a **decision variable**. Continuity of time and of
$(x_1,x_2)$ is imposed at the junction. There is no $u(t)$ to
chatter: each phase is a smooth linear ODE.

Analytic (RGR $=1$, $T=5$, $X_{10}=1$):

$$
t_s^\star = T-1 = 4, \qquad x_2(T)=e^{4}.
$$

The optimizer is **not** given $t_s=4$. Default `tf_guess=5` places
the knot at the midpoint; IPOPT moves it to 4.

---

## 2. `eco_opt` API

```python
from eco_opt import Phase

sys_v = eco_opt.System(states=2, controls=0, model=vegetative_model, name="cohen_veg")
sys_r = eco_opt.System(states=2, controls=0, model=reproductive_model, name="cohen_repro")

p1 = Phase(sys_v, constraints=[(0, "==", t_0(sys_v.t)), ...], control_intervals=20)
p2 = Phase(sys_r, constraints=[(T, "==", t_f(sys_r.t)), ...], control_intervals=20)

ocp = eco_opt.Ocp()
ocp.max(t_f(sys_r.y.repro))
ocp.st(phases=[p1, p2], continuity=True)  # Tf_1 == T0_2 and x continuous
sol = ocp.solve(tf_guess=5.0)
# same transcription, compiled then solved:
# ocp.build(tf_guess=5.0); sol = ocp.optimize()
```

Each phase needs its **own** `System` (distinct CasADi symbols).
`continuity=True` (default) links $t_f^{(k)}=t_0^{(k+1)}$ and
$x^{(k)}(t_k)=x^{(k+1)}(t_k)$ for $\min(n_x^{(k)}, n_x^{(k+1)})$
states.

Results:

- `sol.phases[i]` — `PhaseTrajectory` with `.time`, `.x`, `.u`, `.t0`, `.tf`
- `sol.switch_times` — junction times (`[t_s]` here)
- `sol.t0`, `sol.tf`, `sol.objective` — overall problem
- `sol.time`, `sol.x` — concatenated when every phase has the same $n_x$

`tf_guess` may be a scalar overall $T$ (equal knots), a sequence of
`n_phases+1` knot times, or `n_phases` durations.

---

## 3. What the test checks

`tests/test_multiphase.py`

1. `Phase` is public.
2. Two phases, $t_s$ within $10^{-5}$ of 4, $x_2(T)$ within
   $10^{-6}$ of $e^4$, junction continuity, $x_2\equiv 0$ on
   phase 1, $x_1$ frozen on phase 2.
3. Same numbers from a deliberately wrong knot at $t=2$.

---

## 4. Files

| File | Role |
|---|---|
| `cohen_two_phase.py` | two `Phase`s, baked-in bangs, free $t_s$ |
| `plot_cohen_two_phase.py` | writes `figs/cohen_two_phase.png` |
| `tests/test_multiphase.py` | analytic gate |

---

## 5. Running

```bash
python cohen_two_phase.py
python plot_cohen_two_phase.py
python -m pytest ../../tests/test_multiphase.py -q
```

---

## 6. Limits of this first transcription

- Sequential phases only (no branching).
- No algebraic states, no per-phase parameters.
- Same collocation degree in every phase; interval count may differ.
- `continuity=True` links states by index. Different state layouts
  need `continuity=False` and (not yet) user-supplied link equalities.

A three-phase Goddard (boost / singular / coast) is the natural next
example once custom links are exposed.
