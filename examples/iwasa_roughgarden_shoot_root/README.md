# Iwasa & Roughgarden (1984) shoot / root / fruit allocation

A three-organ annual-plant problem whose analytic optimum is
**bang–singular–bang in two controls**. It is the natural next step
after Cohen (one vegetative compartment, bang-bang) and King &
Roughgarden (one vegetative compartment, mixed-growth singular arc):
here the vegetative body is split into **shoot** and **root**, and the
plant must also choose *when* to send photosynthate to fruit.

The qualitative theorem is Iwasa & Roughgarden (*Theor. Popul. Biol.*
**25**:78–105, 1984). The production function implemented here is the
**Liebig / minimum** special case $G=\min(\alpha S, \beta R)$, for
which the balanced path, the balanced allocations, and the fruit-switch
time are closed form. The three-phase *structure* is theirs; the
specific $G$ is the tractable instance, not a verbatim transcription
of their most general concave production class.

This page is written as example documentation and can be lifted into
the package docs with light editing.

---

## 1. Biological question

A terrestrial plant produces photosynthate with two vegetative organs
that supply complementary resources:

- **Shoot** $S$: light capture, photosynthesis.
- **Root** $R$: water and mineral uptake.

Fruit $F$ is a pure sink: it does not produce. Over a season of known
length $T$ the plant allocates the current photosynthate rate $G(S,R)$
among shoot, root, and fruit so as to maximise season-end fruit $F(T)$.

Iwasa & Roughgarden’s theorem: the optimum is three phases.

1. **Unbalanced bang.** If the seedling is shoot- or root-heavy, put
   *everything* into the currently limiting organ until the plant hits
   the balanced growth path.
2. **Balanced-growth singular arc.** Shoot and root grow together at
   the ratio that maximises daily photosynthesis for the current
   vegetative biomass. Nothing is sent to fruit yet.
3. **Reproductive bang.** Shoot and root stop simultaneously. All
   subsequent photosynthate goes to fruit.

The last phase is Cohen’s flowering switch, applied to the *balanced*
vegetative body. The first phase is the new phenomenon: the plant
corrects its shape before it invests in mixed vegetative growth, and
it will reconverge to the same path after pruning or an environmental
change (the paper’s namesake “shoot/root balance”).

---

## 2. Optimal-control problem (Liebig instance)

$$
\begin{aligned}
\max_{u_S,u_R}\quad
& F(T) \\
\text{s.t.}\quad
& \dot S = u_S\, G, \qquad
  \dot R = u_R\, G, \qquad
  \dot F = (1-u_S-u_R)\, G, \\
& G(S,R) = \min(\alpha S,\, \beta R), \\
& 0 \le u_S,\quad 0 \le u_R,\quad u_S+u_R \le 1, \\
& S(0)=S_0,\; R(0)=R_0,\; F(0)=0, \quad t\in[0,T].
\end{aligned}
$$

The two controls live on a simplex: the unused fraction
$u_F=1-u_S-u_R$ is fruit allocation. There is no separate fruit
control variable.

`eco_opt` transcription: three states, two controls, a *path*
constraint on the derived expression `sys.u[0] + sys.u[1]`, Mayer
objective `ocp.max(t_f(sys.y.F))`.

The true Liebig function `G = min(αS, βR)` is not differentiable on
the balanced path — which is exactly the singular arc. A raw
`cs.fmin` therefore stalls IPOPT (Hessian undefined on a
positive-length interval). The dynamics use a hyperbolic smooth-min

```python
G = 0.5 * (a + b - cs.sqrt((a - b) ** 2 + eps ** 2))  # a=αS, b=βR
dx = cs.vertcat(u_S * G, u_R * G, (1.0 - u_S - u_R) * G)
```

with `eps = 1e-3`. This is a standard C¹ surrogate of `min`, not a
different biological model: tests compare the OCP against the *true*
Liebig closed form. See §8.

---

## 3. Why Liebig, and what is *not* claimed

Iwasa & Roughgarden work with a general production function $G(S,R)$
that is increasing and concave in each argument, with complementary
inputs (more shoot helps only if there is enough root, and vice versa).
Along the balanced path the plant chooses the root/shoot ratio that
maximises $G$ for a given total vegetative biomass
$V=S+R$. That ratio, and the associated allocations, depend on the
shape of $G$.

Liebig $G=\min(\alpha S, \beta R)$ is the extreme complementary case:

- daily photosynthesis is limited by whichever organ is scarcer,
- the maximiser of $G$ at fixed $V$ is the unique ratio that
  *equates* the two arguments, $\alpha S=\beta R$,
- everything after that is algebra.

This is the standard closed-form special case used in teaching notes
and in later shoot/root-allocation papers. It is **not** “the Iwasa &
Roughgarden model with their original $G$ typed in”. The README, the
module docstring, and the figure caption all say so. If a later example
implements a smooth $G$ (e.g. a product or CES form), it should live
beside this one rather than replace it: the Liebig case is the one
with an analytic fruit yield to test against.

---

## 4. Closed-form Liebig solution

Write $\gamma = \alpha\beta/(\alpha+\beta)$ for the balanced relative
growth rate, and

$$
u_S^\star = \frac{\beta}{\alpha+\beta}, \qquad
u_R^\star = \frac{\alpha}{\alpha+\beta}
$$

for the balanced allocations (the unique interior point of the simplex
that keeps $\alpha S=\beta R$ while using all photosynthate
vegetatively).

### 4.1 Balanced growth, once the ratio is right

On $\alpha S=\beta R$ one has $G=\alpha S=\beta R=\gamma V$ with
$V=S+R$, and

$$
\dot S = \gamma S, \qquad \dot R = \gamma R.
$$

Both organs grow exponentially at rate $\gamma$, fruit stays at 0.

### 4.2 Fruit switch (Cohen’s rule with RGR $\gamma$)

After the switch, $S$ and $R$ freeze, so $G$ is constant and
fruit accumulates linearly. The remaining-time argument that gives
Cohen’s $t^\star=T-1$ for RGR $=1$ here gives

$$
t_s = T - \frac{1}{\gamma}.
$$

Terminal fruit is then $F(T)=G(t_s)\,(T-t_s)=G(t_s)/\gamma$.

### 4.3 Unbalanced opening bang

If the seedling is already balanced ($\alpha S_0=\beta R_0$), phase 1
has length zero.

If it is **root-limited** ($\beta R_0 < \alpha S_0$):

- $u_R=1$, $u_S=0$, $G=\beta R$,
- $S$ frozen at $S_0$, $R(t)=R_0 e^{\beta t}$,
- the organs meet when $\beta R=\alpha S_0$, i.e.

$$
t_\mathrm{bal} = \frac1\beta \ln\Bigl(\frac{\alpha S_0}{\beta R_0}\Bigr).
$$

If it is **shoot-limited**, swap the roles of $(\alpha,S)$ and
$(\beta,R)$.

The three-phase picture requires $t_\mathrm{bal}<t_s$. For a very
unbalanced seedling in a short season the plant may never reach the
balanced path before it is forced to fruit; that regime is not the
demo.

### 4.4 Default numbers

$\alpha=\beta=1$ so $\gamma=1/2$, $u_S^\star=u_R^\star=1/2$,
$T=6$, $S_0=1$, $R_0=0.2$ (root-deficient seedling).

$$
t_\mathrm{bal}=\ln 5 \approx 1.60944, \qquad
t_s=4, \qquad
F(T)=\frac{2 e^{2}}{\sqrt{5}}\approx 6.60948.
$$

Trajectory sketch:

| Phase | Time | $u_S,u_R,u_F$ | States |
|---|---|---|---|
| 1 all-to-root | $[0,\ln 5)$ | $0,1,0$ | $S\equiv 1$, $R=0.2 e^{t}$, $F\equiv 0$ |
| 2 balanced | $[\ln 5,\,4)$ | $1/2,1/2,0$ | $S=R=\exp(\tfrac12(t-\ln 5))$ |
| 3 all-to-fruit | $[4,6]$ | $0,0,1$ | $S,R$ frozen, $F=G(t_s)\,(t-4)$ |

---

## 5. Pontryagin sketch (why three phases)

Hamiltonian, maximisation form, with $u_F=1-u_S-u_R$:

$$
H = G\bigl(\lambda_S u_S + \lambda_R u_R + \lambda_F(1-u_S-u_R)\bigr)
  = G\bigl(\lambda_F + u_S(\lambda_S-\lambda_F) + u_R(\lambda_R-\lambda_F)\bigr).
$$

Because $H$ is affine in $(u_S,u_R)$ on a simplex, the optimum is
at a vertex or on an edge/face (singular). The three vertices are
“all shoot”, “all root”, “all fruit”. The shoot–root edge is the
balanced singular arc (both vegetative controls interior, fruit at 0).
Iwasa & Roughgarden identify the switching surfaces with the *marginal
values* of shoot, root, and a unit of fruit, and show:

- if one vegetative organ has strictly higher marginal value than the
  other and than fruit, bang into that organ (phase 1);
- if the two vegetative marginal values are equal and still larger
  than fruit, stay on the balanced edge (phase 2);
- when those common vegetative marginal values fall through 1 (the
  marginal value of immediate fruit), bang to fruit (phase 3).

For Liebig $G$, the balanced edge is exactly $\alpha S=\beta R$,
and the fruit-crossing happens at remaining time $1/\gamma$. The
example does not re-derive the general-$G$ switching surfaces; it
uses the Liebig closed form as an independent check that collocation
recovers the same three phases.

---

## 6. `eco_opt` transcription

```python
sys = eco_opt.System(states=3, controls=2, model=shoot_root_model, name="iwasa_rg")
ocp = eco_opt.Ocp()
ocp.max(eco_opt.t_f(sys.y.F))
ocp.st(
    system=sys,
    constraints=[
        (0.0, "<=", sys.u[0], "<=", 1.0),
        (0.0, "<=", sys.u[1], "<=", 1.0),
        (0.0, "<=", sys.u[0] + sys.u[1], "<=", 1.0),  # path, derived
        (0.0, "<=", sys.y.S),
        (0.0, "<=", sys.y.R),
        (0.0, "<=", sys.y.F),
        (0.0, "==", eco_opt.t_0(sys.t)),
        (T, "==", eco_opt.t_f(sys.t)),
        (S0, "==", eco_opt.t_0(sys.y.S)),
        (R0, "==", eco_opt.t_0(sys.y.R)),
        (F0, "==", eco_opt.t_0(sys.y.F)),
    ],
)
```

API points this example is meant to document:

- **Two controls** with a linear path constraint on a *sum*, not a
  box on a single component. `classify_constraint` treats
  `sys.u[0]+sys.u[1]` as a derived expression, hence a general path
  constraint, while each `sys.u[i]` remains a box.
- **Min-type production** transcribed as a hyperbolic smooth-min of
  Liebig `min(αS, βR)` (`eps=1e-3`), collocation Radau, IPOPT
  `hessian_approximation=limited-memory`. A raw `cs.fmin` puts a
  non-differentiable kink on the singular arc and does not converge.
- **Mayer terminal fruit**, not a Lagrange integral of $u_F G$. The
  two are equivalent along feasible trajectories; Mayer is the
  statement of the biological objective.
- Named signals `sys.y.S`, `sys.y.R`, `sys.y.F`, `sys.y.G`, `sys.y.u_F`.

Initial guess: constant balanced allocations $(u_S,u_R)=(1/2,1/2)$
and the true initial state. The solver is *not* given the three-phase
structure.

---

## 7. What the test checks

`tests/test_iwasa_roughgarden_shoot_root.py`

1. Closed-form identities: $\gamma=1/2$, $t_s=4$,
   $t_\mathrm{bal}=\ln 5$, $F(T)=2e^2/\sqrt{5}$, balanced
   controls $1/2,1/2$.
2. OCP success, and terminal fruit within 3% of the true Liebig value
   (the dynamics use a smooth-min with $\varepsilon=10^{-3}$).
3. Detected $t_\mathrm{bal}$ and $t_s$ within $0.25$ of the
   closed form, in the right order.
4. Phase 1: mean $u_R$ high, mean $u_S$ low before $t_\mathrm{bal}$.
5. Phase 2: mean $|R/S - \alpha/\beta|$ small on the interior of
   the balanced window.
6. Phase 3: $S$ and $R$ nearly frozen after $t_s$, fruit still
   increasing.
7. Simplex: $0 \le u_S+u_R \le 1$ on every interval.

With the smooth-min transcription the default mesh recovers
bang–singular–bang controls that sit on the closed form (all-to-root,
then $u_S=u_R=1/2$, then all-to-fruit). Tests still treat **state
ratios, freeze of $S,R$ after $t_s$, and $F(T)$** as the primary
assertions, so a future mesh that chatters on the simplex edge will
not fail the gate.

---

## 8. Numerical notes

- True Liebig $G=\min(\alpha S,\beta R)$ is nondifferentiable on the
  balanced path. Direct collocation with `cs.fmin` therefore stalls
  (the kink *is* the singular arc). The implemented $G$ is the
  hyperbolic smooth-min of that Liebig function,
  $G_\varepsilon=\tfrac12\bigl(a+b-\sqrt{(a-b)^2+\varepsilon^2}\bigr)$
  with $\varepsilon=10^{-3}$. This is a transcription of the same
  production function, not a CES or product substitute. With these
  defaults the OCP recovers $F(T)$ to $<0.1\%$ of the true Liebig
  closed form, and the three phases in the controls are visually
  bang–singular–bang.
- Radau collocation (same choice as `examples/non_diff/`) with 80
  intervals, degree 3, limited-memory Hessian, `tol=1e-7`.
- The fruit switch is a discontinuity in $u$. Detecting it from
  `u_S+u_R` dropping through $0.4$ is the robust test.
- Do not replace this $G_\varepsilon$ with a CES just to “make it
  smooth”: CES is a different biological model and has no closed-form
  $F(T)$. Keep this file as the analytic Liebig gate.

---

## 9. Relation to other package examples

| Example | Organs | Objective | Optimum |
|---|---|---|---|
| `cohen_king_roughgarden/cohen_terminal_yield.py` | veg + repro | $x_2(T)$ | bang–bang, $t^\star=T-1$ |
| `cohen_king_roughgarden/king_singular_allocation.py` | veg + repro | $\int\ln x_2$ | bang–singular–bang in one control |
| `lindh_2016_flowering/` | veg + repro, three $F(V)$ | $R(T)$ | bang–bang; $t_F^\star(P)$ flips with constraint |
| **this example** | shoot + root + fruit | $F(T)$ | bang–singular–bang in *two* controls on a simplex |
| `microbial_growth/` | resource pool + uptake | $\int G(u)$ | smooth (not bang); **not** a plant-allocation problem |
| `root_opt/` | leaf/root C–N with soil microbes | leaf growth integral | numerical, no closed form |

The Ledder–Manzoni teaching OCP is already `microbial_growth`; it
was deliberately not duplicated.

---

## 10. Files

| File | Role |
|---|---|
| `shoot_root.py` | Liebig OCP, closed-form helpers, trajectory reconstruction |
| `plot_shoot_root.py` | writes `figs/iwasa_roughgarden_shoot_root.png` |
| `tests/test_iwasa_roughgarden_shoot_root.py` | analytic three-phase regression |

Public API: `eco_opt.System`, `eco_opt.Ocp`, `eco_opt.t_0`,
`eco_opt.t_f`, `eco_opt.Namespace`. No MATLAB at runtime.

---

## 11. Parameters

| symbol | default | meaning |
|---|---|---|
| $\alpha$ | 1 | shoot productivity (light) |
| $\beta$ | 1 | root productivity (nutrient/water) |
| $T$ | 6 | season length |
| $S_0$ | 1 | initial shoot |
| $R_0$ | 0.2 | initial root (deficient) |
| $F_0$ | 0 | initial fruit |

Derived: $\gamma=1/2$, $t_s=4$, $t_\mathrm{bal}=\ln 5$,
$F(T)=2e^2/\sqrt{5}$.

---

## 12. Running

```bash
python shoot_root.py
python plot_shoot_root.py
python -m pytest ../../tests/test_iwasa_roughgarden_shoot_root.py -q
```

---

## 13. Mapping onto the `eco_opt` API (for later docs)

| Concept | Where it appears |
|---|---|
| multi-control `System` | `states=3, controls=2` |
| named signals | `sys.y.S`, `.R`, `.F`, `.G`, `.u_F` |
| Mayer objective | `ocp.max(t_f(sys.y.F))` |
| box bounds | `0 <= sys.u[i] <= 1`, states $\ge 0$ |
| derived path constraint | `0 <= sys.u[0]+sys.u[1] <= 1` |
| min-type right-hand side | hyperbolic smooth-min of Liebig $\min(\alpha S,\beta R)$ |
| collocation | Radau, degree 3, 80 intervals |
| IPOPT | `limited-memory` Hessian, `max_iter=2000`, `tol=1e-7` |
| diagnostics | `sol.x` ratios and freeze, not pointwise `sol.u` |

---

## 14. Suggested extensions (not implemented)

- Smooth CES $G(S,R)=\bigl((\alpha S)^{-\rho}+(\beta R)^{-\rho}\bigr)^{-1/\rho}$,
  recovering Liebig as $\rho\to\infty$. No closed-form
  $F(T)$; compare numerically to this file as $\rho$ grows.
- Shoot-deficient seedling ($S_0 < R_0$) — the symmetric opening
  bang. The helpers in `shoot_root.py` already cover it; only the
  default demo and the test use the root-deficient IC.
- Time-varying $\alpha(t),\beta(t)$ (light vs. drought): the
  balanced ratio itself becomes a path, which is a main qualitative
  prediction of the 1984 paper.
- Two-resource (C and N) allocation, as in the 2021 bioRxiv extension
  of Iwasa–Roughgarden. That is a different, larger example.

---

## 15. References

- Iwasa, Y. & Roughgarden, J. (1984). Shoot/root balance of plants:
  optimal growth of a system with many vegetative organs. *Theor.
  Popul. Biol.* **25**:78–105. doi:10.1016/0040-5809(84)90007-8
- Cohen, D. (1971). Maximizing final yield when growth is limited by
  time or by limiting resources. *J. Theor. Biol.* **33**:299–307.
  (single vegetative compartment; fruit switch $T-1/r$)
- King, D. & Roughgarden, J. (1982). Graded allocation between
  vegetative and reproductive growth. *Theor. Popul. Biol.*
  **22**:1–16. (singular mixed growth in one vegetative compartment)
- Iwasa, Y. (2000). Dynamic optimization of plant growth. In
  *The Geometry of Ecological Interactions*, Cambridge.
- Baird, A. et al. (2021). Resource allocation in annual plants.
  bioRxiv 10.1101/2021.04.19.440512. Two-resource (C and N) extension
  of Iwasa–Roughgarden; useful as a future example, not this one.
