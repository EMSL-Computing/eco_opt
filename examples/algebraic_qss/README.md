# Quasi-steady microbial biomass

Litter carbon $C$ is slow. Microbial biomass $M$ is fast, so it is
not given its own $\dot M$. At every instant it sits at
**growth = death**. The control is uptake $u$.

$$
\max_u \int_0^T M\,\mathrm{d}t,\qquad
\dot C=-u,\qquad
0=g(C,M,u),\qquad
C(0)=1,\; T=10,\; 0\le u\le u_{\max}.
$$

$M$ is an **algebraic** unknown: it has no initial condition. The
model returns `(dx, dz, y)` with `algebraics=1`; `dz` is the residual
$g$, enforced as $0=g$.

## Linear closure (`linear_qss.py`)

$g=(e_{\max}-ru)\,u\,C/(K+C)-d M$ is linear in $M$, so there is a
unique formula

$$
M=\frac{(e_{\max}-ru)\,u}{d}\frac{C}{K+C}.
$$

You can substitute that into $J$ (plain ODE) or keep $M$ and
$0=g$. Both must give the same $J$, $C(t)$, and $u(t)$.

## Quadratic death (`quadratic_qss.py`)

Death is $d M + c M^2$, so $g=0$ is a quadratic. The physical root
is the plus branch of a square root; the minus root is junk. Here you
do **not** put `sqrt` in $f$: write $0=g$ and $M\ge 0$.

## Running

```bash
cd examples/algebraic_qss
python linear_qss.py
python quadratic_qss.py
python plot_algebraic_qss.py
python -m pytest ../../tests/test_algebraics.py -q
```

Radau collocation, degree 3.
