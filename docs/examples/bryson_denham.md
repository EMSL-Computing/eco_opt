# Bryson–Denham (`bryson_denham`)

**Question.** Minimum-energy double integrator whose trajectory must
not exceed a path bound on position.

## Formulation

**Dynamics.** \(\dot x=v\), \(\dot v=a\).

**Objective.**

\[
\min \int_0^1 \tfrac12 a^2\,\mathrm{d}t.
\]

**Constraints.** \(x(0)=0\), \(v(0)=1\), \(x(1)=0\), \(v(1)=-1\),
and the path inequality \(x(t)\le 1/9\).

**Horizon.** Fixed, \(T=1\).

Closed form: the trajectory kisses \(x=1/9\) on an interior arc and
\(J^*=4\).

## What was done

Path bound as a constraint on the named position signal. The bound can
be retuned with `set_constraints` / `optimize` without rebuilding the
mesh.

## What we found

\(J\) matches \(4\) to a few parts in \(10^3\). Peak position sits on
\(1/9\). Terminal state matches the specified rest-reversal.

![States](../assets/examples/bryson_denham/bryson_denham_states.png)

![Control](../assets/examples/bryson_denham/bryson_denham_control.png)

## Interpretation

A path constraint is a physiological-style limit: the process may not
leave an admissible set at **any** time, not only at the endpoint.
Collocation enforces that at the nodes.

## Run and tests

```bash
cd examples/bryson_denham
python bryson_denham.py
python -m pytest ../../tests/test_bryson_denham.py -q
```

## References

Bryson and Ho, *Applied Optimal Control* (1975). Also shipped with Yop:
[Bryson–Denham](https://www.yoptimization.com/brysonDenham).
The text and code here are ours.

!!! warning "Not validated against published figures"
    Numerics are checked against the closed-form cost \(J^*=4\), not
    against figures from Bryson and Ho.
