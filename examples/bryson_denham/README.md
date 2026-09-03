# Bryson–Denham

Minimum-energy double integrator with a state path constraint
(Bryson & Ho, *Applied Optimal Control*, 1975).

```
minimize   int_0^1 0.5 * a^2 dt
subject to xdot = v, vdot = a
           x(0)=0, v(0)=1, x(1)=0, v(1)=-1
           x(t) <= 1/9
```

Closed-form cost `J*=4` when the trajectory kisses `x=1/9` on an interior
arc.

To retune the path bound without rebuilding the collocation graph:

```python
ocp, sys = build_ocp()
ocp.build(control_intervals=20, polynomial_degree=5, tf_guess=1.0)
sol_a = ocp.optimize(warm_start=False, tf_guess=1.0)
ocp.set_constraints([(1.0 / 5.0, ">=", sys.y.position)])
sol_b = ocp.optimize()
```

## Running

```bash
python bryson_denham.py
python plot_bryson_denham.py
python -m pytest ../../tests/test_bryson_denham.py -q
```
