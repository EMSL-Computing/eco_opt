# Van der Pol oscillator

The CasADi / PROPT “hello world” for nonlinear Lagrange regulation.
There is no closed-form optimum; the example exists to gate a smooth
oscillator with a quadratic running cost, box bounds on $u$, and a
fixed horizon — the instance every collocation tutorial starts from.

This page is written as example documentation and can be lifted into
the package docs with light editing.

---

## 1. Question

The uncontrolled Van der Pol equation

$$
\ddot x - \mu(1-x^2)\dot x + x = 0
$$

is a nonlinear oscillator with a stable limit cycle (negative damping
near the origin, positive damping at large amplitude). Adding a force
$u$ and a quadratic cost on $(x,\dot x,u)$ asks: how should one
push the oscillator to kill the cycle as cheaply as possible over a
fixed window?

---

## 2. Optimal-control problem

CasADi’s direct-collocation demo (Andersson, Gillis, Diehl) writes the
oscillator in the first-order form

$$
\dot x_1 = (1-x_2^2)x_1 - x_2 + u, \qquad
\dot x_2 = x_1,
$$

which is Van der Pol with $\mu=1$ in the coordinates
$(x_1,x_2)=(\dot x, x)$ up to the sign convention on the damping
term. The OCP is

$$
\begin{aligned}
\min_{u(\cdot)}\quad
& \int_0^{T} \bigl(x_1^2 + x_2^2 + u^2\bigr)\,\mathrm{d}t \\
\text{s.t.}\quad
& \text{dynamics above},\quad
  (x_1,x_2)(0)=(0,1), \\
& -1\le u(t)\le 1, \qquad T=10\text{ fixed}.
\end{aligned}
$$

No terminal-state constraint: the running cost alone drives
$(x_1,x_2)$ to the origin (numerically $\lvert x(T)\rvert<10^{-3}$
on the default mesh). The CasADi demo also boxes $x_1\ge -0.25$; that
path constraint is *not* used here — it is inactive at the unconstrained
optimum of this instance, and omitting it keeps the example a pure
box-on-$u$ Lagrange problem.

PROPT’s textbook listing is the same oscillator as a *Mayer* problem
with an auxiliary integral state $x_3$, $T=5$, and
$-0.3\le u\le 1$. The two are equivalent up to horizon and bounds.
This file follows the CasADi collocation demo because `eco_opt` is a
CasADi transcription; the PROPT Mayer form is the same pattern as
Goddard’s integral-via-`t_f(x3)`.

---

## 3. Pontryagin (no closed form)

Hamiltonian:

$$
H = x_1^2+x_2^2+u^2
  + \lambda_1\bigl((1-x_2^2)x_1-x_2+u\bigr)
  + \lambda_2 x_1.
$$

Affine in the *dynamics*, quadratic in the *cost*, so the unconstrained
minimiser is the linear feedback

$$
u^\star = -\tfrac12\lambda_1,
$$

projected onto $[-1,1]$. On this instance the upper bound $u=1$ is
active at the beginning (the oscillator is kicked hard toward the
origin) and the control then rides the interior formula. There is no
singular arc ($H_{uu}=2\neq 0$).

Costates are not returned by `eco_opt`; the tests do not reconstruct
them. The bound-saturation of $u$ is the cheap PMP check.

---

## 4. `eco_opt` transcription

```python
sys = eco_opt.System(states=2, controls=1, model=vdp_model, name="vdp")
ocp = eco_opt.Ocp()
ocp.min([sys.y.x1 ** 2 + sys.y.x2 ** 2 + sys.y.u ** 2])
ocp.st(
    system=sys,
    constraints=[
        (-1.0, "<=", sys.u[0], "<=", 1.0),
        (0.0, "==", eco_opt.t_0(sys.t)),
        (10.0, "==", eco_opt.t_f(sys.t)),
        (0.0, "==", eco_opt.t_0(sys.y.x1)),
        (1.0, "==", eco_opt.t_0(sys.y.x2)),
    ],
)
```

A Python list in `ocp.min([...])` is a running (Lagrange) cost. Mayer
would wrap `t_f(...)`. Collocation: 40 Legendre intervals, degree 4,
constant guess $u=0$.

Open-loop replay (`ocp.simulate(sol)`) integrates $\dot x=f(t,x,u^*)$
with CVODES and should track `sol.x`. Shrinking-horizon
`ocp.closed_loop(x0, dt, t_final=T)` on the **same** $f$ reprints
that tape. A plant with an extra term in $\dot x_1$ does not.

---

## 5. What the test checks

There is no analytic $J^\star$. The default-mesh objective is
**regression-tracked**:

    REGRESSION_OBJECTIVE = 2.8878976747752554

computed with `control_intervals=40`, `polynomial_degree=4`, Legendre.
A silent change in the transcription that moves $J$ by more than
$10^{-6}$ fails the test.

`tests/test_van_der_pol.py`

1. Solver success, $T=10$, ICs.
2. $J$ matches the frozen regression value.
3. $u$ stays in $[-1,1]$; $\max u > 0.95$ (upper bound active).
4. $x(T)\approx 0$ (the cost regulates without a terminal constraint).
5. A coarser mesh ($N=20$, degree 3) agrees on $J$ to 3%.

Mesh note: $J$ still drops a few $10^{-2}$ from $N=20$ to $N=80$
because the initial bang $u=1$ is a discontinuity. The regression
number is tied to the default mesh on purpose; do not “improve” it
without updating the constant.

---

## 6. Files

| File | Role |
|---|---|
| `van_der_pol.py` | OCP, frozen regression objective |
| `plot_van_der_pol.py` | writes `figs/van_der_pol.png` |
| `tests/test_van_der_pol.py` | regression + bound + mesh check |

Public API: `eco_opt.System`, `eco_opt.Ocp`, `eco_opt.t_0`,
`eco_opt.t_f`, `eco_opt.Namespace`. No MATLAB at runtime.

---

## 7. Parameters

| symbol | default | meaning |
|---|---|---|
| $T$ | 10 | fixed horizon |
| $(x_1,x_2)(0)$ | $(0,1)$ | CasADi demo IC |
| $u$ | $[-1,1]$ | additive force |

---

## 8. Running

```bash
cd examples/van_der_pol
python van_der_pol.py
python plot_van_der_pol.py
python -m pytest ../../tests/test_van_der_pol.py -q
```

---

## 9. Mapping onto the `eco_opt` API (for later docs)

| Concept | Where it appears |
|---|---|
| Lagrange objective | `ocp.min([x1**2 + x2**2 + u**2])` |
| fixed horizon | `t_f(sys.t) == 10` |
| box on control | `-1 <= sys.u[0] <= 1` |
| free terminal state | no `t_f` state equalities |
| named signals | `sys.y.x1`, `.x2`, `.u` |
| collocation | Legendre, degree 4, 40 intervals |

Related examples: Lenhart 17.3 (LQ tracking, also Lagrange, linear
plant); Xuezhong Wang ex.1 (LQ double integrator, closed form). Van der
Pol is the nonlinear analogue: same cost structure, no analytic $u(t)$.

---

## 10. Suggested extensions (not implemented)

- PROPT Mayer form: auxiliary $x_3$ with $\dot x_3=x_1^2+x_2^2+u^2$,
  $T=5$, $-0.3\le u\le 1$, minimise $x_3(T)$.
- Drive-to-origin as a hard terminal constraint
  $x_1(T)=x_2(T)=0$ (Andersson 2012 paper, slightly different $J$).
- Limit-cycle tracking rather than regulation.

---

## 11. References

- Van der Pol, B. (1926). On relaxation-oscillations.
  *Phil. Mag.* **2**:978–992.
- Andersson, J., Gillis, J., Horn, G., Rawlings, J. B. & Diehl, M.
  (2019). CasADi: a software framework for nonlinear optimization and
  optimal control. *Math. Prog. Comp.* **11**:1–36. Direct collocation
  demo: `docs/examples/python/direct_collocation.py` in the CasADi
  tree (this instance: $T=10$, $x(0)=(0,1)$, $u\in[-1,1]$,
  $L=x_1^2+x_2^2+u^2$).
- PROPT manual, “Van der Pol Oscillator”: Mayer form with $T=5$ and
  tighter control bounds. Same plant, different horizon/bounds.
