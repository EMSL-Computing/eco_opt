"""Open-loop dry-forecast irrigation: S, u, and sinusoidal ET."""
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import irrigation_bucket as ib

FIG_DIR = os.path.join(HERE, "figs")


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    sol, sys, _ = ib.solve()
    t = sol.time
    S = sol.x[:, 0]
    u = np.append(sol.u[:, 0], sol.u[-1, 0])
    E = ib.E0 * (1.0 + ib.E_AMP * np.sin(2.0 * np.pi * t))

    fig, axes = plt.subplots(3, 1, figsize=(6.5, 7.2), sharex=True)
    axes[0].plot(t, S, "k-")
    axes[0].axhline(ib.S_MIN, color="0.5", ls="--", lw=0.8)
    axes[0].set_ylabel("S (mm)")
    axes[0].grid(alpha=0.3)
    axes[1].step(t, u, where="post", color="C1")
    axes[1].set_ylabel("u (mm d$^{-1}$)")
    axes[1].grid(alpha=0.3)
    axes[2].plot(t, E, color="C2")
    axes[2].set_ylabel("E(t) (mm d$^{-1}$)")
    axes[2].set_xlabel("t (d)")
    axes[2].grid(alpha=0.3)
    fig.suptitle(
        f"Dry-forecast irrigation OCP  J=∫u dt={sol.objective:.2f} mm  "
        f"S(T)={S[-1]:.1f}",
        fontsize=10,
    )
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "irrigation_bucket.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
