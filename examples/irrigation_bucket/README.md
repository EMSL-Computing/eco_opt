# Irrigation bucket

A one-state root-zone water balance, inspired by the irrigation MPC of
Shang, Chen, Stroock & You (2020). This is **not** a reproduction of
their robust forecasts or the 40% water-saving number.

The controller believes it is a dry fortnight with known sinusoidal
evapotranspiration. Rain, if any, happens only on the **plant**.

$$
\dot S = u + P(t) - E(t) - k S,
\qquad
\min \int_0^T u\,dt,
\qquad
S_{\min}\le S\le S_{\max}.
$$

$S$ is storage (mm), $u$ irrigation (mm d$^{-1}$).
$E(t)=E_0\bigl(1+a\sin(2\pi t)\bigr)$ with period 1 d (always positive).
The OCP uses $P\equiv 0$.

## Open vs closed loop

On the dry plant, shrinking-horizon closed loop reprints the open-loop
tape (same $f$, no disturbance).

On a plant that rains on days 5–7, the tape keeps watering; closed loop
sees the jump in $S$ and cuts $u$. Water use falls; $S$ stays off
the ceiling.

## Run

```bash
python irrigation_bucket.py
python plot_irrigation_bucket.py
python open_vs_closed.py
```

Figures: `figs/irrigation_bucket.png`, `figs/open_vs_closed.png`.

## Parameters

| symbol | value | meaning |
|---|---|---|
| $T$ | 14 d | fixed horizon |
| $S(0)$ | 40 mm | initial storage |
| $S_{\min},S_{\max}$ | 20, 80 mm | floor / saturation |
| $k$ | 0.05 d$^{-1}$ | linear drainage |
| $E_0,a$ | 4 mm d$^{-1}$, 0.6 | ET mean and daily amplitude |
| $u_{\max}$ | 15 mm d$^{-1}$ | irrigation cap |
| rain | 10 mm d$^{-1}$ on $[5,7]$ | plant only |

## Citation (inspiration only)

Shang, C., Chen, W.-H., Stroock, A. D. & You, F. (2020). Robust model
predictive control of irrigation systems with active uncertainty
learning and data analytics. *IEEE Transactions on Control Systems
Technology* 28: 1493–1504. https://doi.org/10.1109/TCST.2019.2916753
