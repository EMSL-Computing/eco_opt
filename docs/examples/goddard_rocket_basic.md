# Goddard, nondimensional (`goddard_rocket_basic`)

**Question.** Same ascent problem in nondimensional variables.

Transcription check, not ecology.

## Formulation

Maximise final altitude, free \(T\). Constants \(m_0=1\), \(m_f=0.6\),
\(h_0=1\), \(g_0=1\), \(D_0=310\), \(c=0.5\), \(\beta=500\),
\(T_{\max}=3.5\). Distinct from the dimensional archive in
`goddard_benchmark`.

**Objective.** \(\max h(T)\). **Horizon.** Free. Thrust and mass
bounds.

## What was done

Nondimensional `System`; free \(t_f\); compared loosely to literature
BVP values (drag-law offset vs some COPS listings).

## What we found

A bang–singular–coast thrust; altitude slightly above 1; \(T\) near
\(0.2\). Tests use a few-percent tolerance plus mesh refinement.

![States](../assets/examples/goddard_rocket_basic/goddard_rocket_basic_states.png)

![Control](../assets/examples/goddard_rocket_basic/goddard_rocket_basic_control.png)

## Interpretation

A second scaling of Goddard is useful as a mesh gate. It is not a
second ecological application.

## Run and tests

```bash
cd examples/goddard_rocket_basic
python goddard_rocket_basic.py
python -m pytest ../../tests/test_goddard_rocket_basic.py -q
```

## References

Related Yop page:
[Basic rocket](https://www.yoptimization.com/goddardRocket).
The text and code here are ours.
