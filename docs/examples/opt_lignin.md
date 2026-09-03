# Ligninolytic investment (`opt_lignin`)

**Question.** Over a fixed decomposition period, when should a
decomposer community pay the carbon cost of oxidising aromatic litter?

## Formulation

| Symbol | Role | Meaning |
|---|---|---|
| \(C_T(t)\) | state | total remaining litter C |
| \(C_O(t)\) | state | aromatic (lignin-derived) C |
| \(v_O(t)\) | control | ligninolytic oxidation rate |
| \(L=C_O/C_T\) | diagnostic | aromatic fraction |

**Dynamics.**

\[
\begin{aligned}
D_H &= v_{H,\max}\,g(L)\,(C_T-C_O),
\qquad
D_O = v_O C_O,\\
e &= e_{\max}-r_O v_O,
\qquad
G = e\,(D_H+D_O),\\
\dot C_T &= -D_H-D_O+G,
\qquad
\dot C_O &= -D_O+\mu_O G,\\
g(L)&=\exp\bigl(-(L/a)^b\bigr).
\end{aligned}
\]

Growth efficiency \(e\) falls linearly with ligninolytic effort.

**Objective.**

\[
\max_{v_O(\cdot)}\; J=\int_0^{T} G(t)\,\mathrm{d}t.
\]

**Constraints.** \(C_T(0)\), \(C_O(0)\) given; \(C_T\ge 0\), \(C_O\ge 0\);
\(0\le v_O\le v_{O,\max}\). Horizon \(T\) **fixed**. Kinetic coefficients
\(v_{H,\max}\), \(\mu_O\), \(r_O\) are pinned parameters.

## What was done

The published process ODE is reused. \(v_O(t)\) is the only control.
Parameters are `System` parameters retuned with `set_constraints` /
`optimize` on a frozen mesh (sweep of cost \(r_O\) and related
coefficients).

## What we found

Optimal \(v_O\) is near zero, then rises after a lag and later declines:
hydrolytic substrate is used first. Higher \(r_O\) delays onset and
lowers the peak. Baseline (\(r_O=25\,\mathrm{d}\), litter C:N 50)
gives \(J^*=0.1628\).

![Litter C and ligninolytic rate](../assets/examples/opt_lignin/opt_lignin_figure2.png)

## Interpretation

A time-invariant oxidation rate hides the lag. The hypothesis “maximise
growth net of oxidative cost over a fixed incubation” is enough to
produce delayed ligninolysis. That lag is a prediction of the
hypothesis, not an extra switch put into \(f\).

## Open loop vs closed loop

`open_vs_closed.py` replays \(v_O^*\) on the same litter ODE, then on a
process that receives extra litter C in days 600–900. On the forecast
model the shrinking closed loop reprints the tape. After the unforecast
C pulse the tape keeps the old ligninolytic schedule; the closed loop
re-solves from the measured \(C_T\) and can raise \(v_O\). Same lesson
as microbial uptake and irrigation: the hypothesis is unchanged; the
**state** from which it is solved is not.

![Open vs closed loop](../assets/examples/opt_lignin/open_vs_closed.png)

## Run and tests

```bash
cd examples/opt_lignin
python lignin_decomposition.py
python open_vs_closed.py
python -m pytest ../../tests/test_opt_lignin.py -q
```

## References

Chakrawal et al. (2024),
[doi:10.1111/nph.19572](https://doi.org/10.1111/nph.19572).
