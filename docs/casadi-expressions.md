# CasADi expressions

The model callback is traced **once** when `System` is built. Arguments
`t`, `x`, `u`, and (if used) `z` and `p` are CasADi `MX` symbols, not
Python floats or NumPy arrays. IPOPT later evaluates that graph.
`math.exp` and `np.exp` do not work on those symbols.

After `solve`, `sol.time`, `sol.x`, `sol.u` **are** NumPy. Plotting
and post-processing use NumPy and matplotlib
([Results and plots](results-and-plots.md)).

```python
import casadi as cs
import eco_opt
from eco_opt import Namespace

k = 0.3


def model(t, x, u):
    stock, harvest = x[0], u[0]
    # stock is MX: cs.exp, not np.exp or math.exp
    loss = k * stock * cs.exp(-stock)
    y = Namespace(stock=stock, harvest=harvest, loss=loss)
    return -harvest - loss, y
```

Python `+`, `-`, `*`, `/`, `**` on `MX` mixed with Python floats are
valid. Return a Python list of derivatives; `System` stacks it to
match `states`. A one-state model may return a scalar.

Numeric bounds in constraint tuples stay Python floats
(`0`, `1.0 / 9.0`, `float("inf")`). Expressions in the model, in
`ocp.max` / `ocp.min`, and in constraint *expressions* are CasADi.

A Python `if` on a **numeric** constant (a float you pass in, not
`sys.p`) is evaluated before the graph exists. A Python `if` on
`t`, `x`, `u`, `z`, or `p` is invalid
([Transcription](transcription.md#what-may-appear-in-f)).

Import `casadi as cs` only when the right-hand side needs a function
below. Many one-state models need only `+`, `-`, `*`, `/`.

## Arithmetic

| | |
|---|---|
| Add, subtract, multiply, divide | `a + b`, `a - b`, `a * b`, `a / b` |
| Power | `a ** b` or `cs.power(a, b)` |
| Square | `a ** 2` |
| Unary minus | `-a` |

Do not stack \(\dot{\mathbf{x}}\) with `cs.vertcat`. Return
`[dx1, dx2, ...]`.

## Elementary

| | |
|---|---|
| Exponential | `cs.exp(x)` |
| Natural log | `cs.log(x)` |
| Base-10 log | `cs.log10(x)` |
| Square root | `cs.sqrt(x)` |
| Absolute value | `cs.fabs(x)` |
| Sign (\(-1,0,+1\)) | `cs.sign(x)` |
| Floor, ceiling | `cs.floor(x)`, `cs.ceil(x)` |
| Remainder | `cs.fmod(x, y)` |

`cs.exp` and `cs.sign` appear in Goddard; `cs.log` in King allocation,
shallow-lake \(J\), and a Lenhart lab; `cs.sqrt` in the Iwasa smooth
minimum and quadratic QSS.

Floor, ceiling, and `fmod` are stepwise. Prefer a smooth surrogate if
IPOPT’s Hessian must exist (see below).

## Trigonometric and hyperbolic

| | |
|---|---|
| Sine, cosine, tangent | `cs.sin`, `cs.cos`, `cs.tan` |
| Inverse | `cs.asin`, `cs.acos`, `cs.atan`, `cs.atan2(y, x)` |
| Hyperbolic | `cs.sinh`, `cs.cosh`, `cs.tanh` |
| Inverse hyperbolic | `cs.asinh`, `cs.acosh`, `cs.atanh` |
| \(\pi\) | `cs.pi` |

`cs.sin` / `cs.cos` appear in the brachistochrone, isoperimetric
example, and irrigation ET. `cs.tanh` is a useful smooth step even
though no shipped example uses it.

## Min, max, error function

| | |
|---|---|
| Minimum, maximum | `cs.fmin(a, b)`, `cs.fmax(a, b)` |
| Error function | `cs.erf(x)` |
| Inverse error function | `cs.erfinv(x)` |

`cs.fmin` / `cs.fmax` are nondifferentiable on a set of positive
measure. Lindh uses `cs.fmax(V, 1e-9)` as a guard; lignin uses
`cs.fmin` on a **numeric** efficiency cap (evaluated before the
graph). If IPOPT stalls on a `min` / `max` **of states**, replace it
with a smooth surrogate. [Shoot–root](examples/iwasa_roughgarden_shoot_root.md)
uses

\[
\tfrac12\bigl(a+b-\sqrt{(a-b)^2+\varepsilon^2}\bigr)
\]

in place of \(\min(a,b)\). A smooth step is `0.5 * (1 + cs.tanh(x / eps))`.

## Branching and logic

A switch that depends on \(t\), \(\mathbf{x}\), or \(\mathbf{u}\) is
`cs.if_else`, not a Python `if`.

| | |
|---|---|
| Switch | `cs.if_else(cond, a, b)` |
| And / or / not | `cs.logic_and`, `cs.logic_or`, `cs.logic_not` |
| Comparisons | `t >= t0`, `x[0] > 0` (overloaded; use inside `if_else`) |
| \(\pm\infty\) in **tuples** | `float("inf")` (or `cs.inf`) |

```python
alive = cs.if_else(t >= t_cut, 0.0, 1.0)
rate = alive * vmax * x[0] / (km + x[0])

on = cs.logic_and(t >= t0, t < tf)
I = I0 * cs.if_else(on, 1.0, 0.0)
```

Pulse forcing: [non_diff](examples/non_diff.md), irrigation,
microbial and lignin `open_vs_closed`. Nested `if_else` in
[root_opt](examples/root_opt.md) selects mineralisation vs
immobilisation.

## Tabulated forcing

Weather, light, or a prescribed input that is a known function of
time is **not** a control. Store it as a CasADi interpolant of \(t\):

```python
I_fun = cs.interpolant("I_light", "linear", [t_grid], I_tab)
I = I_fun(t)
```

See [greenhouse](examples/greenhouse.md). Do not call `np.interp`
inside the model.

## Linear algebra (rarely needed)

Process ODEs are almost always scalar or a short list of
scalars. If a term is genuinely a matrix expression, CasADi also
provides `cs.mtimes`, `cs.inv`, `cs.solve`, `cs.dot`, `cs.norm_2`,
`cs.sum1`, `cs.transpose`. Prefer a Python list of scalar rates
over building one large `MX` column yourself.

## Reference

CasADi’s own list of operations:
[web.casadi.org/docs](https://web.casadi.org/docs/).
