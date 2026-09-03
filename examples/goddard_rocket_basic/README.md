# Goddard rocket (nondimensional, max ascent)

COPS-style nondimensional Goddard: maximize final altitude, free final
time. Distinct from the dimensional free-tf problem in
`examples/goddard_benchmark/`.

Parameters: `m0=1`, `mf=0.6`, `h0=1`, `g0=1`, `D0=310`, `c=0.5`,
`beta=500`, `Tmax=3.5`.

Literature (Petit et al., IFAC 2008, indirect BVP): `h(T)=1.01283692`,
`T=0.19885626`. The drag law here omits the `h0` offset used in some
COPS listings, so the test uses a few-percent tolerance plus mesh
refinement.

## Running

```bash
python goddard_rocket_basic.py
python plot_goddard_rocket_basic.py
python -m pytest ../../tests/test_goddard_rocket_basic.py -q
```
