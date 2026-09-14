# Greenhouse climate control

van Straten, van Ooteghem, van Willigenburg & van Henten, *Optimal
Control of Greenhouse Cultivation*. Crop dry weight and indoor
temperature; heating is the control. Light and outside temperature are
**exogenous**: they enter $f$ as functions of $t$, not as controls.
They are stored as linear interpolants of a two-day table (clipped
daylight sinusoid). That is the point of the example: tabulated
weather in an OCP, with no extra NLP variables.

## Problem

$$
\begin{aligned}
\min_u\quad
& -p_5\,x_1(T)+\int_0^{T} p_4 u\,\mathrm{d}t \\
\dot x_1 &= p_1 I(t)\,x_2, \qquad
\dot x_2 = p_2\bigl(T_{\mathrm{out}}(t)-x_2\bigr)+p_3 u, \\
x_1(0)&=0,\; x_2(0)=10,\; T=48\,\mathrm{h},\; 0\le u\le 10.
\end{aligned}
$$

$x_1$ is crop dry weight [kg m$^{-2}$], $x_2$ indoor temperature
[°C], $u$ heating [W m$^{-2}$]. $I(t)$ is incoming light (zero at
night), $T_{\mathrm{out}}(t)$ outside temperature. Minimising that
$J$ maximises crop value minus heating cost.

No closed form. Default-mesh $J$ is regression-tracked.

`greenhouse.py` builds the interpolants from an in-code table.
`greenhouse_from_file.py` reads the same series from `forcing.csv`
(`t,light,T_out`). The OCP is identical; the test checks that $J$,
$x$, and $u$ match.

## Running

```bash
python greenhouse.py
python greenhouse_from_file.py
python plot_greenhouse.py
python -m pytest ../../tests/test_greenhouse.py -q
```
