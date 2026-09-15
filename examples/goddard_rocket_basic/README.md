# Goddard rocket (nondimensional, max ascent)

COPS-style nondimensional Goddard: maximize final altitude, free final
time. Distinct from the dimensional free-tf problem in
`examples/goddard_benchmark/`.

Parameters: `m0=1`, `mf=0.6`, `h0=1`, `g0=1`, `D0=310`, `c=0.5`,
`beta=500`, `Tmax=3.5`.

Literature (Petit et al., IFAC 2008, indirect BVP): `h(T)=1.01283692`.
Height is the comparison (few-percent tolerance plus mesh refinement).
The drag law here is `D0*exp(-beta*h)` (no `h0` offset), unlike
`D0*exp(-beta*(h-h0)/h0)` in some COPS/Petit listings. This
transcription's `T` is not the indirect-BVP `T=0.19885626` (drag-law /
piecewise-constant `u` / mesh).

## Running

```bash
cd examples/goddard_rocket_basic
python goddard_rocket_basic.py
python plot_goddard_rocket_basic.py
python -m pytest ../../tests/test_goddard_rocket_basic.py -q
```
