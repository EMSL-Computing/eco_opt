# Nondifferentiable pulse forcing

Fixed-horizon Mayer problem with a discontinuous `if_else` pulse in the
dynamics. No `y` signal tree. Radau collocation, degree 3.

```
minimize   x3(tf)
subject to x1dot = x2
           x2dot = -x1 - x2 + u + d(t)
           x3dot = 5*x1^2 + 2.5*x2^2 + 0.5*u^2
           x(0)=0, tf=2
```

`d(t)` is a rectangular pulse on `[0.5, 0.6]`. No closed form; the
objective is regression-tracked.

## Running

```bash
python non_diff.py
python plot_non_diff.py
python -m pytest ../../tests/test_non_diff.py -q
```
