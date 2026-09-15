# First solve

The dimensional Goddard rocket: maximise final altitude with free
burnout time. This is a **traditional optimal-control problem**. The
script below is complete. Ecological models use the same structure;
they are listed under [Examples](examples/index.md).

Full write-up, figures, and the archived reference:
[Goddard, dimensional](examples/goddard_benchmark.md).

If IPOPT fails, see [Troubleshooting](troubleshooting.md). For what
counts as a state, a control, or a cost, see
[Posing a problem](posing.md).
The control is piecewise constant on the collocation mesh; see
[Transcription](transcription.md).
This problem has a **free** horizon (`tf_guess=200`); time bounds are
written explicitly in `st`
([Constraints and time](constraints-and-time.md)).

```python
import casadi as cs
import eco_opt


def goddard_model(t, x, u):
    v, h, m = x[0], x[1], x[2]
    F = u[0]

    D0 = 0.01227
    beta = 0.145e-3
    c = 2060.0
    g0 = 9.81
    r0 = 6.371e6

    D = D0 * cs.exp(-beta * h)
    F_D = cs.sign(v) * D * v ** 2
    g = g0 * (r0 / (r0 + h)) ** 2

    dv = (F * c - F_D) / m - g
    dh = v
    dm = -F
    dx = cs.vertcat(dv, dh, dm)

    y = eco_opt.Namespace()
    y.rocket = eco_opt.Namespace(
        velocity=v, height=h, mass=m, fuelMassFlow=F
    )
    y.drag = eco_opt.Namespace(coefficient=D, force=F_D)
    y.gravity = g
    return dx, y


sys = eco_opt.System(states=3, controls=1, model=goddard_model, name="goddard")
rocket = sys.y.rocket

m0 = 214.839
mf = 67.9833
Fm = 9.525515

ocp = eco_opt.Ocp()
ocp.max(eco_opt.t_f(rocket.height))
ocp.st(
    system=sys,
    constraints=[
        (0, "==", eco_opt.t_0(sys.t)),
        (0, "==", eco_opt.t_0(rocket.velocity)),
        (0, "==", eco_opt.t_0(rocket.height)),
        (m0, "==", eco_opt.t_0(rocket.mass)),
        (0, "<=", eco_opt.t_f(sys.t), "<=", float("inf")),
        (0, "<=", rocket.velocity, "<=", float("inf")),
        (0, "<=", rocket.height, "<=", float("inf")),
        (mf, "<=", rocket.mass, "<=", m0),
        (0, "<=", rocket.fuelMassFlow, "<=", Fm),
    ],
)
sol = ocp.solve(
    control_intervals=59,
    polynomial_degree=5,
    collocation_points="legendre",
    tf_guess=200.0,
)
print(sol.success, sol.tf, sol.objective)
```

Same solve as a module: `examples/goddard_benchmark/goddard_benchmark.py`.

`sol.success` should be true. This mesh gives \(T\approx 213\,\mathrm{s}\)
and \(h(T)\approx 1.61\times 10^{5}\,\mathrm{m}\). The example page
compares the same transcription to an independent archive.
