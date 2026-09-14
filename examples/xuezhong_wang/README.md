# Double-integrator examples

Two problems from Wang, X. (n.d.). Solving optimal control problems
with MATLAB — Indirect methods, 21.

## ex.1 — minimum-energy double integrator, fixed final time

    minimize   int_0^2 0.5*a^2 dt
    subject to xdot = v, vdot = a
               x(0)=0, v(0)=0, x(2)=5, v(2)=2

Closed form: `a*(t) = -4.5*t + 5.5`, `J* = 7.75`.

## ex.2 — nonnegative control, fixed final time

    minimize   int_0^2 0.5*a^2 dt
    subject to xdot = v, vdot = a, a >= 0
               x(0)=1, v(0)=1, x(2)=3, v(2)>=0

Trivial optimum: `u*(t)=0`, `J*=0`.

## Running

```bash
python bryson_denham_fixed_tf.py
python double_integrator_nonneg_control.py
python plot_xuezhong_wang.py
pytest ../../tests/test_xuezhong_wang.py -q
```
