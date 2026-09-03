# Cohen (1971) and King & Roughgarden (1982) plant allocation

Two closely related annual-plant optimal-control problems on the *same*
two-compartment growth model, differing only in the fitness criterion.
Together they are the cleanest plant-allocation pair with closed-form
controls that `eco_opt` currently ships:

| Criterion | Paper | Control structure | Closed form |
|---|---|---|---|
| Terminal seed yield \(x_2(T)\) | Cohen (1971) | bang–bang | switch at \(T-1\) |
| \(\int_0^T \ln x_2(t)\,dt\) | King & Roughgarden (1982) | bang–singular–bang | singular ratio \(x_1=(T-1-t)x_2\) |

They sit next to `examples/lindh_2016_flowering/` (Cohen’s model with
three production functions, no log-running-payoff) and
`examples/iwasa_roughgarden_shoot_root/` (two vegetative organs plus
fruit). They are *not* the same problem as
`examples/microbial_growth/` (Ledder–Manzoni substrate uptake).

This page is written as example documentation: biological question,
optimal-control statement, Pontryagin derivation, how the problem is
posed in `eco_opt`, what the tests actually check, and numerical
caveats. It can be lifted into the package docs with light editing.

---

## 1. Biological question

An annual plant has a finite season \([0, T]\) and partitions
photosynthate between a **vegetative body** \(x_1\) (leaves, stems,
roots that produce more photosynthate) and a **reproductive body**
\(x_2\) (flowers, fruits, seeds that do not). Instantaneous production
is taken to be linear in vegetative mass, with relative growth rate
scaled to 1:

\[
\dot x_1 = u\, x_1, \qquad
\dot x_2 = (1-u)\, x_1, \qquad
0 \le u(t) \le 1.
\]

Here \(u\) is the fraction of photosynthate *kept* in the vegetative
body. The remainder \(1-u\) is sent to reproduction.

Two fitness criteria have been used for this model:

1. **Cohen (1971).** In a season of known length, selection maximises
   final seed yield \(x_2(T)\). The optimum is an instantaneous switch
   from pure vegetative growth to pure reproduction — the classical
   “bang-bang flowering time”.
2. **King & Roughgarden (1982).** In a fluctuating environment the
   relevant long-run fitness is the geometric mean of yield, which in
   continuous time is represented by a running payoff \(\ln x_2\).
   Concavity of \(\ln\) makes mixed vegetative/reproductive growth
   (a *singular* arc) optimal for part of the season: a graded
   allocation, not a single switch.

Cohen’s result is why so many phenology models are bang-bang. King &
Roughgarden’s result is why mixed growth is not a modelling bug: it is
the PMP-optimal response to a concave-in-reproduction payoff.

---

## 2. Shared plant model

States: vegetative mass \(x_1\), reproductive mass \(x_2\).
Control: vegetative allocation \(u \in [0,1]\).
Time: fixed horizon \([0, T]\).

\[
\dot x_1 = u\, x_1, \qquad
\dot x_2 = (1-u)\, x_1.
\]

The second equation is **not** \(\dot x_2 = (1-u)x_2\). Photosynthate
is produced by the vegetative body. Atkins, Aghaee, Martcheva & Hager
(arXiv:2010.06744) reproduce this problem as their Example 2; the HTML
rendering of their eq. (70) writes \(x_2\) on the right-hand side, but
their problem statement (72) and the rest of the derivation use \(x_1\).
This example follows the consistent (Cohen / King–Roughgarden / (72))
form.

In `eco_opt` the model is a two-state, one-control `System`:

```python
def plant_model(t, x, u):
    veg, repro = x[0], x[1]
    alloc = u[0]
    dx = cs.vertcat(alloc * veg, (1.0 - alloc) * veg)
    y = Namespace(veg=veg, repro=repro, alloc=alloc)
    return dx, y
```

The two examples differ only in the objective and in the initial
reproductive mass (Cohen allows \(x_2(0)=0\); the log-payoff requires
\(x_2(0)>0\)).

---

## 3. Cohen (1971): Mayer terminal yield, bang-bang

### 3.1 Problem

\[
\max_{u(\cdot)}\; x_2(T)
\quad\text{s.t.}\quad
\dot x_1 = u x_1,\;
\dot x_2 = (1-u) x_1,\;
x_1(0)=X_{10}>0,\;
x_2(0)=0,\;
0\le u\le 1.
\]

Defaults: \(T=5\), \(X_{10}=1\).

### 3.2 Pontryagin derivation

Hamiltonian (maximisation form):

\[
H = \lambda_1 u x_1 + \lambda_2 (1-u) x_1
  = \lambda_2 x_1 + u x_1 (\lambda_1 - \lambda_2).
\]

Switching function \(\sigma = x_1(\lambda_1-\lambda_2)\). Since
\(x_1>0\),

\[
u^\star =
\begin{cases}
1 & \lambda_1 > \lambda_2,\\
0 & \lambda_1 < \lambda_2,\\
\text{singular} & \lambda_1 = \lambda_2.
\end{cases}
\]

Costates: \(\dot\lambda_1 = -u(\lambda_1-\lambda_2)-\lambda_2\),
\(\dot\lambda_2 = 0\). Free terminal states give
\(\lambda_1(T)=\lambda_2(T)=0\), so \(\lambda_2\equiv 0\). Then
\(\dot\lambda_1 = -u\lambda_1\) with \(\lambda_1(T)=0\), which forces
\(\lambda_1\equiv 0\) *after* the last time \(u=1\). The remaining-time
calculation with RGR \(=1\) puts the unique switch at

\[
t^\star = T-1.
\]

(If the relative growth rate were \(r\) rather than 1, the switch would
be at \(T-1/r\). That is the same rule Lindh et al. (2016) use for
exponential production, and the same rule the Iwasa–Roughgarden Liebig
example uses with balanced RGR \(\gamma\).)

Integrating the bang-bang dynamics:

\[
\begin{aligned}
x_1(t) &= X_{10}e^{t},
  && 0\le t\le t^\star,\\
x_1(t) &= X_{10}e^{t^\star},
  && t^\star\le t\le T,\\
x_2(t) &= 0,
  && 0\le t\le t^\star,\\
x_2(t) &= X_{10}e^{t^\star}(t-t^\star),
  && t^\star\le t\le T,\\
x_2(T) &= X_{10}e^{T-1}.
\end{aligned}
\]

With equal conversion efficiencies, \(x_1(T)=x_2(T)\): the plant
switches when remaining time is just enough for the frozen vegetative
capital to produce its own mass in seed.

Defaults \(T=5\), \(X_{10}=1\) give \(t^\star=4\) and
\(x_2(T)=e^4 \approx 54.598\).

There is no singular arc: \(\lambda_1=\lambda_2=0\) on a positive-length
interval would contradict the maximisation of \(x_2(T)\) (one can always
move a little allocation earlier). Direct collocation, given a
continuous relaxation \(u\in[0,1]\), saturates at the bounds with a
single switch — the same numerical test used in
`examples/lindh_2016_flowering/` for Theorem 1 of Lindh et al.

### 3.3 `eco_opt` transcription

Mayer objective on the named reproductive signal:

```python
ocp = eco_opt.Ocp()
ocp.max(eco_opt.t_f(sys.y.repro))
ocp.st(
    system=sys,
    constraints=[
        (0.0, "<=", sys.u[0], "<=", 1.0),
        (0.0, "==", eco_opt.t_0(sys.t)),
        (T, "==", eco_opt.t_f(sys.t)),
        (x10, "==", eco_opt.t_0(sys.y.veg)),
        (0.0, "==", eco_opt.t_0(sys.y.repro)),
        ...
    ],
)
```

Collocation: 80 Legendre intervals, polynomial degree 3. Constant
initial guess \(u=0.5\) — the solver is *not* told that the control is
bang-bang.

### 3.4 What the test checks

`tests/test_cohen_king_roughgarden.py::test_cohen_matches_analytic_bang_bang`

- solver success
- \(x_2(T)\) within \(5\times 10^{-2}\) of \(e^{4}\)
- flowering time (last interval with \(u>0.5\)) within \(0.15\) of \(4\)
- \(u(0)\approx 1\), \(u(T)\approx 0\)
- \(x_1(T)\approx x_2(T)\) (equal-efficiency corollary)

A solver reporting success is not enough; the numbers have to match the
closed form.

---

## 4. King & Roughgarden (1982): Lagrange log-yield, singular arc

### 4.1 Problem

\[
\max_{u(\cdot)}\; \int_0^T \ln x_2(t)\,dt
\quad\text{s.t.}\quad
\dot x_1 = u x_1,\;
\dot x_2 = (1-u) x_1,\;
x_1(0)=X_{10}>0,\;
x_2(0)=X_{20}>0,\;
0\le u\le 1.
\]

The running payoff \(\ln x_2\) is the continuous analogue of Cohen’s
discrete geometric-mean criterion (King & Roughgarden, *Theor. Popul.
Biol.* **22**:1–16, 1982; Atkins et al. 2020, Example 2). Because
\(\ln x_2\) is undefined at 0, the demo starts with \(X_{20}=1>0\).

Defaults used here: \(T=5\), \(X_{20}=1\), and the plant is started
*on the singular manifold*

\[
X_{10} = (T-1)X_{20} = 4,
\]

so the initial vegetative bang is absent (\(t_1=0\)). That is the
cleanest analytic gate: one singular arc joined to a terminal
reproductive bang.

### 4.2 Pontryagin derivation

Hamiltonian (now a maximisation of \(\int \ln x_2\)):

\[
H = \ln x_2 + \lambda_1 u x_1 + \lambda_2 (1-u) x_1
  = \ln x_2 + \lambda_2 x_1 + u x_1(\lambda_1-\lambda_2).
\]

Switching function \(\sigma = x_1(\lambda_1-\lambda_2)\). Costates:

\[
\dot\lambda_1 = -u(\lambda_1-\lambda_2)-\lambda_2, \qquad
\dot\lambda_2 = -1/x_2,
\]

with free-endpoint transversality \(\lambda_1(T)=\lambda_2(T)=0\).

On a singular interval \(\lambda_1=\lambda_2=\lambda\), so
\(\dot\lambda_1=-\lambda\) and \(\dot\lambda_2=-1/x_2\), hence
\(\lambda=1/x_2\). Differentiating \(\lambda=1/x_2\) and matching
\(\dot\lambda=-\lambda\) gives \(\dot x_2 = x_2\). Combined with the
state equation \(\dot x_2=(1-u)x_1\),

\[
u^\star(t) = 1 - \frac{x_2}{x_1}.
\]

The ratio \(\rho=x_1/x_2\) then satisfies \(\dot\rho=-1\), so \(\rho(t)\)
is affine. Matching the terminal reproductive bang (next paragraph)
fixes the intercept:

\[
x_1(t) = (T-1-t)\, x_2(t), \qquad
u^\star(t) = 1 - \frac{1}{T-1-t}
\quad\text{on }[t_1,t_2).
\]

The singular arc cannot continue all the way to \(T\) (the factor
\(T-1-t\) would change sign, and transversality \(\lambda_2(T)=0\) is
incompatible with \(\lambda=1/x_2\neq 0\)). It joins a terminal bang
\(u\equiv 0\) on \([t_2,T]\). On that bang, \(x_1\) is frozen and \(x_2\)
grows linearly. Continuity of \(\lambda_1,\lambda_2\) and of their first
derivatives at \(t_2\) produces a *parameter-free* nonlinear system in

\[
y = T-t_2, \qquad z = x_2(t_2)/x_1(t_2).
\]

The system reduces to a scalar transcendental equation for
\(w=1/z\):

\[
e^w = w^2 + w + 1, \qquad w>0
\]

(\(w=0\) is the trivial root). Then

\[
y = 1+w, \qquad z = 1/w, \qquad
\frac{x_2(T)}{x_1(T)} = y+z.
\]

Newton on the residual gives the published values (Atkins et al.,
eqs. 96–99)

\[
y \approx 2.79328213, \quad
z \approx 0.55763674, \quad
\frac{x_2(T)}{x_1(T)} \approx 3.35091887.
\]

Hence \(t_2 = T-y \approx 2.2067\) for \(T=5\), independent of the
initial masses. The initial masses only decide whether an opening
vegetative bang is needed to *reach* the manifold.

If \(X_{10}/X_{20} > T-1\), phase 1 is \(u=1\) until
\(x_1(t_1)=(T-1-t_1)x_2(t_1)\). That variant is implemented analytically
in the derivation but is not the default demo (starting on the manifold
is the stricter numeric test of the singular arc itself).

### 4.3 Closed-form states, start-on-manifold

On \([0,t_2]\):

\[
x_2(t) = X_{20} e^{t}, \qquad
x_1(t) = (T-1-t) X_{20} e^{t}, \qquad
u(t) = 1 - \tfrac{1}{T-1-t}.
\]

On \([t_2,T]\): \(u=0\), \(x_1\) frozen, \(x_2\) linear. The running
payoff integrates in closed form:

\[
J^\star = T\ln X_{20} + \tfrac12 t_2^2 + y\, t_2 + w^2.
\]

With \(X_{20}=1\) the log term drops. This is the number the OCP
objective is compared against.

### 4.4 `eco_opt` transcription

Lagrange objective — a Python list means “integrate this expression”:

```python
ocp.max([cs.log(sys.y.repro)])
```

Reproductive mass is boxed away from zero (`x2 >= 1e-8`) so that
\(\ln x_2\) stays defined. Collocation is denser than Cohen’s (120
Legendre intervals, degree 4) because the control is singular.

### 4.5 What the test checks — and what it deliberately does not

Direct collocation on a singular arc typically **chatters**: \(u(t)\)
oscillates around the analytic singular value instead of sitting on it.
That is a known artefact of discretising an affine-in-control Hamiltonian
(Atkins et al. 2020 discuss it at length; TV regularisation is one
remedy, not used here). The tests therefore treat

- the **state ratio** \(x_1/((T-1-t)x_2)\approx 1\) on the interior of
  \([0,t_2]\) (this is the singular-arc test; after \(t_2\) the same
  ratio leaves 1 and later returns through 1 as \(T-1-t\to 0\), so a
  last-time-near-1 detector is the wrong junction finder),
- the **terminal ratio** \(x_2(T)/x_1(T)\approx 3.3509\),
- the **objective** \(J\),
- the junction \(t_2\) from the first drop of \(u\) to the reproductive
  bang,
- and a late-time reproductive bang (\(u\approx 0\) after \(t_2\))

as the diagnostics. Pointwise \(u(t)\) on the singular window is plotted
against the analytic curve; with the default mesh it actually sits on
the closed form, but the tests do not require that.

`tests/test_cohen_king_roughgarden.py` also checks that the Newton
junction constants reproduce the published \(y,z\) to \(10^{-8}\).

---

## 5. How the two problems differ, in one picture

Cohen’s payoff only cares about \(x_2\) at the *end* of the season, so
there is no reason to grow reproductive tissue early: every unit of
vegetative capital compounds until the last moment it can still pay for
itself. King & Roughgarden’s payoff scores \(\ln x_2\) *every day*, so
an early seed set has value, but so does keeping some vegetative capital
to grow future \(x_2\). The compromise is the mixed-growth singular arc,
followed by the same “all remaining photosynthate to seed” bang that
Cohen uses.

If you only need a bang-bang flowering-time demo, Cohen (or Lindh
exponential) is enough. If you need a plant problem whose analytic
optimum is *not* bang-bang, King–Roughgarden is the standard example.

---

## 6. Files

| File | Role |
|---|---|
| `cohen_terminal_yield.py` | Cohen Mayer OCP, closed-form switch/yield/trajectories |
| `king_singular_allocation.py` | King Lagrange OCP, junction Newton, singular-ratio helpers |
| `plot_cohen_king.py` | writes `figs/cohen_terminal_yield.png` and `figs/king_singular_allocation.png` |
| `tests/test_cohen_king_roughgarden.py` | analytic regression for both |

Public API used: `eco_opt.System`, `eco_opt.Ocp`, `eco_opt.t_0`,
`eco_opt.t_f`, `eco_opt.Namespace`. No MATLAB at runtime.

---

## 7. Parameters

**Cohen.** \(T=5\), \(X_{10}=1\), \(x_2(0)=0\), \(0\le u\le 1\).
Analytic: \(t^\star=4\), \(x_2(T)=e^4\).

**King.** \(T=5\), \(X_{20}=1\), \(X_{10}=4\) (on-manifold), \(0\le u\le 1\),
\(x_2\ge 10^{-8}\). Analytic: \(t_1=0\), \(t_2=T-y\approx 2.2067\),
\(x_2(T)/x_1(T)\approx 3.3509\).

---

## 8. Running

From this directory, with the package installed:

```bash
python cohen_terminal_yield.py
python king_singular_allocation.py
python plot_cohen_king.py
python -m pytest ../../tests/test_cohen_king_roughgarden.py -q
```

---

## 9. Numerical notes

- Cohen is easy: affine control, unique switch, Mayer objective. Default
  collocation is already close to machine-precision on \(x_2(T)\).
- King is the first package example whose *analytic* optimum is
  singular. Expect chatter in \(u\) on \([0,t_2)\). Refining the mesh
  reduces the state-ratio error; it does not necessarily produce a
  visually constant singular control. That is a property of the
  transcription, not a failure of IPOPT.
- Do not “fix” chatter with an example-specific penalty on \(\dot u\)
  unless the same regularisation is exposed as a general option. The
  tests are written to survive chatter.
- A constant initial guess \(u=0.6\) is enough when the plant starts on
  the manifold. Off-manifold ICs may need a better guess or more
  intervals.
- Collocation points: Legendre, matching the other analytic plant
  examples (Lindh, Lenhart 17.2).

---

## 10. Mapping onto the `eco_opt` API (for later docs)

| Concept | Cohen | King |
|---|---|---|
| `System(states, controls, model)` | 2, 1, `plant_model` | same |
| named signals | `sys.y.veg`, `sys.y.repro` | same |
| objective | Mayer: `ocp.max(t_f(sys.y.repro))` | Lagrange: `ocp.max([cs.log(sys.y.repro)])` |
| control bounds | box on `sys.u[0]` | same |
| time horizon | `t_0(sys.t)==0`, `t_f(sys.t)==T` | same |
| initial state | `t_0` equalities | same, with \(x_2(0)>0\) |
| path bound | \(x_1,x_2\ge 0\) | \(x_2\ge 10^{-8}\) for \(\ln\) |
| result | `sol.objective`, `sol.x`, `sol.u`, `sol.time` | plus state-ratio diagnostic |

The list-wrapped objective is the package convention for a running cost
(see `examples/xuezhong_wang/bryson_denham_fixed_tf.py`,
`examples/lenhart_workman/linear_bangbang_ex17_2.py`). A `t_f(...)`
wrapped objective is Mayer (see Goddard, Lindh, Cohen).

---

## 11. Suggested extensions (not implemented)

- Off-manifold King ICs: \(X_{10}/X_{20} > T-1\), with \(t_1\) from
  \(X_{10}e^{t_1}=(T-1-t_1)X_{20}\).
- Relative growth rate \(r\neq 1\): replace \(T-1\) by \(T-1/r\).
- King & Roughgarden (1982, *TPB* **21**:194–204), the companion paper
  on *multiple* bang switches in a deterministic environment with
  intermediate seed set — different problem, same plant.
- TV-regularised singular control (Atkins et al.), as a package-level
  option rather than an example hack.

---

## 12. References

- Cohen, D. (1971). Maximizing final yield when growth is limited by
  time or by limiting resources. *J. Theor. Biol.* **33**:299–307.
- King, D. & Roughgarden, J. (1982). Graded allocation between
  vegetative and reproductive growth for annual plants in growing
  seasons of random length. *Theor. Popul. Biol.* **22**:1–16.
  doi:10.1016/0040-5809(82)90032-6
- King, D. & Roughgarden, J. (1982). Multiple switches between
  vegetative and reproductive growth in annual plants. *Theor. Popul.
  Biol.* **21**:194–204. (companion paper; not this example)
- Atkins, S., Aghaee, M., Martcheva, M. & Hager, W. (2020). Solving
  singular control problems in mathematical biology, using PASA.
  arXiv:2010.06744. Example 2 is the King–Roughgarden plant; junction
  constants \(y,z\) used here come from their eqs. 96–99.
- Lindh, M. et al. (2016). Constrained growth flips the direction of
  optimal phenological responses among annual plants. *New Phytol.*
  **209**:1591–1599. Same Cohen plant with three \(F(V)\); see
  `examples/lindh_2016_flowering/`.
- Iwasa, Y. & Roughgarden, J. (1984). Shoot/root balance of plants.
  *Theor. Popul. Biol.* **25**:78–105. Two vegetative organs; see
  `examples/iwasa_roughgarden_shoot_root/`.
