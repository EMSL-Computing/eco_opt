# Isoperimetric constraint

Lagrange problem with an auxiliary integral state enforcing
`int_0^1 u^2 dt = 10`. No `y` signal tree.

```
minimize   int_0^1 x1 dt
subject to x1dot = -sin(x1) + u,  x2dot = u^2
           x1(0)=1, x2(0)=0, x1(1)=0, x2(1)=10
           -4 <= u <= 4
```

No published reference for this instance; the objective is
regression-tracked.

## Running

```bash
python iso_constraint.py
python plot_iso_constraint.py
python -m pytest ../../tests/test_iso_constraint.py -q
```
