# Goddard rocket (dimensional, max ascent)

Free-final-time maximum-altitude Goddard in SI units. Same vehicle as
`examples/multiphase_goddard_landing/` (that example adds a landing
phase).

```
maximize   h(tf)
subject to vdot = (F*c - sign(v)*D(h)*v^2)/m - g(h)
           hdot = v,  mdot = -F
           (v,h,m)(0)=(0,0,m0),  mf <= m <= m0,  0 <= F <= Fm
```

`m0=214.839`, `mf=67.9833`, `Fm=9.525515`. Compared against the archived
trajectory in `work/goddard_benchmark/yop_goddard.mat`.

## Running

```bash
cd examples/goddard_benchmark
python goddard_benchmark.py
python plot_goddard_benchmark.py
python -m pytest ../../tests/test_goddard_benchmark.py -q
```
