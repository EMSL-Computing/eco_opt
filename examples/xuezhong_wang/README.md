# xuezhong_wang

Python/`eco_opt` port of two examples from Wang, X. *Solving optimal
control problems with MATLAB — Indirect methods* (source dir:
`examples/OptLignin/YOP benchmark/Xuezhong Wang_Solving optimal control
problems with MATLAB/`).

## ex.1 — minimum-energy double integrator, fixed final time
(`bryson_denham_fixed_tf.py`)

    minimize   int_0^2 0.5*a^2 dt
    subject to xdot = v, vdot = a
               x(0)=0, v(0)=0, x(2)=5, v(2)=2

Port of `Xuezhong_Wang_ex_1.m` (25 control intervals in the original).
This is a linear-quadratic two-point boundary value problem with a known
closed-form solution: the costate is affine in t, so the optimal control
is affine, `a*(t) = -4.5*t + 5.5`, giving analytic objective `J* = 7.75`
(full derivation in the module docstring).

Validation: `tests/test_xuezhong_wang.py::test_ex1_matches_analytic_lq_solution`
checks the eco_opt objective and boundary control values against the
closed form (rel. tol ~5e-3 at 100 control intervals, degree 5).

## ex.2 — double integrator with nonnegative control, fixed final time
(`double_integrator_nonneg_control.py`)

    minimize   int_0^2 0.5*a^2 dt
    subject to xdot = v, vdot = a, a >= 0
               x(0)=1, v(0)=1, x(2)=3, v(2)>=0

Port of `Xuezhong_Wang_ex_2.m` (50 control intervals in the original).
The `a>=0` sign constraint plus the free (inequality-only) terminal
velocity makes the *trivial* constant-velocity coast (`u*(t)=0`,
`v(t)=1`) exactly feasible and optimal: it satisfies `x(2)=1+1*2=3` and
`v(2)=1>=0` with zero control effort, so `J*=0`.

Validation: `tests/test_xuezhong_wang.py::test_ex2_trivial_zero_control_optimum`
checks `J≈0`, `u(t)>=0` everywhere, and the boundary conditions.

## Running

```bash
python bryson_denham_fixed_tf.py
python double_integrator_nonneg_control.py
python plot_xuezhong_wang.py   # writes figs/xuezhong_wang_ex{1,2}.png
pytest ../../tests/test_xuezhong_wang.py -q
```

Figures compare the eco_opt control/state trajectories against the
analytic reference (ex.1) or annotate the trivial optimum (ex.2); see
`figs/`.
