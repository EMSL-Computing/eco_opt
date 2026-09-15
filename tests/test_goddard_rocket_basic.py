import os
import sys as _sys

HERE = os.path.dirname(__file__)
_sys.path.insert(0, os.path.join(HERE, "..", "examples", "goddard_rocket_basic"))

from goddard_rocket_basic import solve, REFERENCE_HEIGHT, MF


def test_goddard_rocket_basic_physical_sanity_and_convergence():
    sol, sys = solve(control_intervals=100)
    assert sol.success

    # Physical sanity: fuel fully spent to mf, monotonic-ish height increase
    # to a peak, thrust saturates near Tmax then drops to 0.
    assert abs(sol.x[-1, 2] - MF) < 1e-3
    assert sol.x[:, 1].max() > sol.x[0, 1]  # height increased
    assert sol.tf > 0

    # Reference literature value (Petit et al. IFAC 2008 indirect BVP shoot,
    # same nondimensional Goddard parameter set): h(T) = 1.01283692.
    # Drag law here is D0*exp(-beta*h) (no h0 offset), which differs
    # slightly from the COPS/Petit convention D0*exp(-beta*(h-h0)/h0);
    # this explains the ~1-2% gap from the literature height. We use a
    # generous tolerance and rely primarily on grid-refinement below.
    # Do not assert sol.tf against the BVP T=0.19885626: piecewise-constant
    # u, mesh, and the drag-law offset make this transcription's T different.
    rel_err = abs(sol.objective - REFERENCE_HEIGHT) / REFERENCE_HEIGHT
    assert rel_err < 0.03, f"objective={sol.objective} vs reference {REFERENCE_HEIGHT}"

    # Grid-refinement convergence: objective should stabilize as N grows.
    sol50, _ = solve(control_intervals=50)
    sol200, _ = solve(control_intervals=200)
    assert sol50.success and sol200.success
    rel_change = abs(sol200.objective - sol.objective) / abs(sol.objective)
    assert rel_change < 1e-3, f"objective not converged: N=100 {sol.objective}, N=200 {sol200.objective}"
