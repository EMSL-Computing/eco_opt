# Shallow-lake eutrophication

A minimal model of a **clear vs turbid lake** and of how much phosphorus
loading agriculture should be allowed to add. Internal recycling from
sediments is a positive feedback: once phosphorus is high, the lake
recycles still more phosphorus. That can produce **alternative stable
states** and **hysteresis** (Scheffer et al. 1993; Carpenter, Ludwig
& Brock 1999). The optimal-control problem is a different object: for
some parameters there are two competing long-run management paths
(Mäler, Xepapadeas & de Zeeuw 2003; a Skiba-type split in the sense of
Skiba 1978 and Wagener 2003). Those two ideas are not the same.

    ecological tipping point  ≠  Skiba / two NLP locals

The first lives in $\dot P=f(P,L)$ at constant loading. The second
lives in the whole OCP (objective + discounting + $P(0)$).

---

## Ecology

Phosphorus in the water column $P$, external loading $L$ (control):

$$
\dot P
=
\underbrace{L}_{\text{external loading}}
-
\underbrace{b P}_{\text{export / sedimentation}}
+
\underbrace{\frac{P^{2}}{1+P^{2}}}_{\text{internal recycling}}.
$$

This is the dimensionless recycling used by Mäler et al. (2003), the
Hill form $r P^{q}/(m^{q}+P^{q})$ with $q=2$, $m=1$, $r=1$
(Carpenter et al. 1999). When $P$ is small, recycling is weak. When
$P$ is large, recycling saturates near 1 and keeps the lake dirty.

For $b=0.52$ and a constant loading in a window
($L\approx 0.02$–$0.07$) there are **three** equilibria:

$$
P_{\mathrm{clear}}\ \text{stable},
\qquad
P_{\mathrm{threshold}}\ \text{unstable},
\qquad
P_{\mathrm{turbid}}\ \text{stable}.
$$

At $L=0.05$: $P\approx 0.13$, $0.58$, $1.31$.

Panel (a) of the figure is that hysteresis curve (loading vs $P$).
Raising $L$ past the upper fold collapses a clear lake to turbid.
Lowering $L$ back to the collapse value does **not** restore it: you
must go below the lower fold (hysteresis).

---

## Where optimal control enters

Loading is not just “bad”. It stands for economically useful activity
(agriculture). The canonical planner problem of Mäler, Xepapadeas &
de Zeeuw (2003) is, on a long finite horizon,

$$
\max_{L(\cdot)}
\int_0^{T} e^{-\rho t}
\bigl(
\underbrace{\ln L}_{\text{benefit of loading}}
-
\underbrace{c P^{2}}_{\text{damage from a dirty lake}}
\bigr)\,\mathrm{d}t.
$$

Question: how much loading to allow, given that too much can lock the
lake in a persistent turbid state?

Defaults here: $b=0.52$, $\rho=0.03$, $c=1$, $T=80$, $P(0)=1$,
$10^{-3}\le L\le 1.5$.

**Two constant guesses for the same $P(0)$** converge to different
locals (both `Solve_Succeeded`):

| $L_{\mathrm{guess}}$ | $J$ (max) | $P(T)$ |
|---|---|---|
| $0.05$ | $\approx -135.22$ | $\approx 2.00$ (dirtier) |
| $0.50$ | $\approx -132.85$ | $\approx 1.54$ (cleaner, better $J$) |

That is panel (b). It is **not** a time-varying $L(t)$ seed: a scalar
high vs low loading is enough. It is the reason this example exists —
initialization can pick which long-run lake you manage toward.

A true Skiba point (Skiba 1978; Wagener 2003) in the infinite-horizon
theory is an initial $P$ where the two values are *equal*. Here the
two $J$s differ, so these are two locals, one better.

---

## Restoration (second OCP)

Start turbid ($P(0)=1.8$) and require $P(T)\le P_{\mathrm{clean}}=0.4$,
while staying close to a desired long-run loading $L_{\mathrm{desired}}=0.05$
(inside the hysteresis window, on the clear branch):

$$
\min_{L(\cdot)}\ \int_0^{T} \bigl(L-L_{\mathrm{desired}}\bigr)^{2}\,\mathrm{d}t.
$$

Because of hysteresis, holding $L\equiv L_{\mathrm{desired}}$ from a
turbid start need not reach $P_{\mathrm{clean}}$. The solver cuts
loading hard (across the restore fold), then can ease toward
$L_{\mathrm{desired}}$. Panel (c). On this mesh a constant guess
already finds that path; the ecology is still the tipping point, not
the seed.

---

## Files

| File | Role |
|---|---|
| `shallow_lake.py` | ODE, both OCPs, equilibria |
| `plot_shallow_lake.py` | writes `figs/shallow_lake.png` |
| `tests/test_shallow_lake.py` | three equilibria; management split; restoration terminal |

## Running

```bash
python shallow_lake.py
python plot_shallow_lake.py
python -m pytest ../../tests/test_shallow_lake.py -q
```

## References

- Scheffer, M., Hosper, S. H., Meijer, M.-L., Moss, B. & Jeppesen, E.
  (1993). Alternative equilibria in shallow lakes. *Trends in Ecology
  & Evolution* **8**:275–279.
  https://doi.org/10.1016/0169-5347(93)90254-M
- Carpenter, S. R., Ludwig, D. & Brock, W. A. (1999). Management of
  eutrophication for lakes subject to potentially irreversible change.
  *Ecological Applications* **9**:751–771.
  https://doi.org/10.1890/1051-0761(1999)009[0751:MOEFLS]2.0.CO;2
- Mäler, K.-G., Xepapadeas, A. & de Zeeuw, A. (2003). The economics of
  shallow lakes. *Environmental and Resource Economics* **26**:603–624.
  https://doi.org/10.1023/B:EARE.0000007351.99227.42
- Wagener, F. O. O. (2003). Skiba points and heteroclinic
  bifurcations, with applications to the shallow lake system. *Journal
  of Economic Dynamics and Control* **27**:1533–1561.
  https://doi.org/10.1016/S0165-1889(02)00070-2
- Skiba, A. K. (1978). Optimal growth with a convex-concave production
  function. *Econometrica* **46**:527–539.
  https://doi.org/10.2307/1914227

