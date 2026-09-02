# Lenhart & Workman examples

Three problems from Lenhart & Workman, *Optimal Control Applied to
Biological Models*.

## ex.17.2 — linear objective, bang-bang control

    maximize   int_0^2 (2*x - 3*u) dt
    subject to xdot = x + u, 0 <= u <= 2
               x(0) = 5

Analytic switch time `t_switch = 2 - ln(2.5) ~= 1.0837`.

## ex.17.3 — least-squares tracking, free terminal state

    minimize   int_0^2 (x - t^2)^2 dt
    subject to xdot = u, 0 <= u <= 4
               x(0) = 1

Smooth LQ tracking (not bang-bang). Regression-tracked objective.

## ex.20.4 — minimum-time trolley stop, free terminal time

    minimize   tf
    subject to xdot = v, vdot = a, -1 <= a <= 1
               x(0) = 10, v(0) = 0, x(tf) = 0, v(tf) >= 0

Analytic `tf* = 2*sqrt(x0) ~= 6.32456` for `x0=10`.

## Running

```bash
python linear_bangbang_ex17_2.py
python lq_tracking_ex17_3.py
python min_time_trolley_ex20_4.py
python plot_lenhart_workman.py
pytest ../../tests/test_lenhart_workman.py -q
```
