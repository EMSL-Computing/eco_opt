# Standard result plots

Use these conventions when visualization will help a user understand an
eco_opt solution. Plotting is optional and uses matplotlib from the
`[examples]` extra; it is not part of eco_opt's core API.

## Default figure

Produce a compact stack of time-series panels:

1. one panel per state, drawn against `sol.time`;
2. one panel per control, drawn as piecewise-constant stairs using
   `sol.time_intervals`; and
3. derived outputs, algebraics, constraint boundaries, or comparison
   trajectories only when they help answer the user's question.

Always use scientific names and units in axis labels. Share the time axis,
place its label on the bottom panel, use a light grid, and avoid legends when
only one unnamed trajectory is present. Include the objective and final time
in the title when they are meaningful.

Do not plot controls as ordinary lines through `sol.time[:-1]`. Controls are
piecewise constant, so use `Axes.stairs(values, sol.time_intervals)`.

Use this adaptable baseline:

```python
from pathlib import Path

import matplotlib.pyplot as plt


def plot_solution(sol, state_labels, control_labels, output=None):
    """Plot an OcpResults trajectory using eco_opt grid semantics."""
    if not sol.success:
        raise RuntimeError(sol.return_status)

    nx = sol.x.shape[1]
    nu = sol.u.shape[1] if sol.u.ndim == 2 else 0
    if len(state_labels) != nx or len(control_labels) != nu:
        raise ValueError("provide one label per state and control")

    nrows = nx + nu
    fig, axes = plt.subplots(
        nrows,
        1,
        figsize=(7.0, max(3.0, 2.2 * nrows)),
        sharex=True,
        squeeze=False,
    )
    axes = axes[:, 0]

    for i, label in enumerate(state_labels):
        axes[i].plot(sol.time, sol.x[:, i], color="C0", linewidth=2)
        axes[i].set_ylabel(label)

    for j, label in enumerate(control_labels):
        ax = axes[nx + j]
        ax.stairs(sol.u[:, j], sol.time_intervals, color="C1", linewidth=2)
        ax.set_ylabel(label)

    for ax in axes:
        ax.grid(alpha=0.25)

    axes[-1].set_xlabel("time")
    fig.suptitle(f"Optimal solution: J={sol.objective:.4g}, tf={sol.tf:.4g}")
    fig.tight_layout()

    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=150, bbox_inches="tight")

    return fig, axes
```

Adapt rather than blindly copy the template. For a model with no states or no
controls, omit the corresponding panels and construct a nonempty figure from
the signals that are present.

## Useful additions

- Plot a named or derived signal with `sol.signal(expression)` against
  `sol.time`.
- Plot algebraic variables from `sol.z` against `sol.time` and label them as
  algebraic quantities.
- Overlay `ocp.simulate(sol)` with a dashed line when replay accuracy or mesh
  error matters; label the collocation and simulation trajectories clearly.
- Mark state or control bounds with thin horizontal dashed lines when active
  constraints are scientifically important.
- For multiple phases, either plot each `sol.phases[i]` explicitly or mark
  `sol.switch_times`; do not hide discontinuities by connecting incompatible
  phase variables.
- For closed-loop comparisons, use consistent colors for the same quantity
  and different line styles for open loop, replay, and closed loop.

## Delivery

In an interactive session, return the figure and let the caller decide when to
show it. In a noninteractive script, select the `Agg` backend before importing
`matplotlib.pyplot`, save the figure, and close it after saving. Use a clear
filename and report its location.

Treat plots as diagnostics as well as presentation. Point out visible bound
activity, switching, discontinuities, replay disagreement, or suspicious
oscillation, but do not infer biological validity from appearance alone.
