# lenhart_workman

Python/`eco_opt` port of three examples from Lenhart, S. & Workman, J. T.
*Optimal Control Applied to Biological Models* (source dir:
`examples/OptLignin/YOP benchmark/Examples from Lenhart and
Workman_Optimal control applied to biological models/`).

## ex.17.2 — linear objective, bang-bang control
(`linear_bangbang_ex17_2.py`, port of `BangBangControl_lenhart_ex17_2.m`)

    maximize   int_0^2 (2*x - 3*u) dt
    subject to xdot = x + u, 0 <= u <= 2
               x(0) = 5

The Hamiltonian is affine in `u`, so Pontryagin's maximum principle gives
bang-bang control with a single switch. Solving the costate ODE
analytically (module docstring) gives the exact switching time
`t_switch = 2 - ln(2.5) ~= 1.0837`.

Validation: `tests/test_lenhart_workman.py::test_ex17_2_bangbang_switch_matches_analytic`
checks the eco_opt switch time against this closed form (tol 0.1, at
30 control intervals) and that `u` saturates at the bounds on each side.

## ex.17.3 — least-squares tracking, free terminal state
(`lq_tracking_ex17_3.py`, port of `BangBangControl_lenhart_ex17_3.m` —
despite the filename this is a smooth tracking problem, not bang-bang)

    minimize   int_0^2 (x - t^2)^2 dt
    subject to xdot = u, 0 <= u <= 4
               x(0) = 1

No simple closed form; validated via a regression-tracked objective and
sanity checks that `u` stays within its box bounds while `x(t)` tracks
the reference parabola `t^2` (see figure).

## ex.20.4 — minimum-time trolley stop, free terminal time
(`min_time_trolley_ex20_4.py`, port of
`Free Terminal Time Problems_ex20_4.m`)

    minimize   tf
    subject to xdot = v, vdot = a, -1 <= a <= 1
               x(0) = 10, v(0) = 0, x(tf) = 0, v(tf) >= 0

Classic minimum-time double-integrator "drive to rest" problem: bang-bang
with one switch (full deceleration toward the origin, then full
acceleration to arrive at rest exactly at x=0). By symmetry the analytic
minimum time is `tf* = 2*sqrt(x0) ~= 6.32456` for `x0=10` — this matches
the MATLAB script's own in-code comment
(`disp("analytical terminal time is " + num2str(2*x0^0.5))`).

Validation: `tests/test_lenhart_workman.py::test_ex20_4_min_time_matches_analytic`
checks `tf` against the closed form (tol 1e-2) and confirms the bang-bang
control shape (`u(0)<-0.9`, `u(tf)>0.9`).

## Running

```bash
python linear_bangbang_ex17_2.py
python lq_tracking_ex17_3.py
python min_time_trolley_ex20_4.py
python plot_lenhart_workman.py   # writes figs/lenhart_workman_ex{17_2,17_3,20_4}.png
pytest ../../tests/test_lenhart_workman.py -q
```
