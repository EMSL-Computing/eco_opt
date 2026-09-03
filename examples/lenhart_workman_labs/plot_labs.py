"""Figures for Lenhart and Workman labs."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))
import lab1_introductory as lab1
import lab2_mold_fungicide as lab2
import lab3_bacteria as lab3
import lab4_bounded as lab4
import lab5_cancer as lab5
import lab6_fish as lab6
import lab7_epidemic as lab7
import lab8_hiv as lab8
import lab9_bears as lab9
import lab10_glucose as lab10
import lab11_timber as lab11
import lab12_bioreactor as lab12
import lab13_predator_prey as lab13

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def _save(fig, name):
    fig.tight_layout()
    out = os.path.join(FIGDIR, name)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("wrote", out)


def plot_lab1():
    sol, _ = lab1.solve()
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0])
    ax[0].set_ylabel("$x$")
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[1].set_ylabel("$u$")
    ax[1].set_xlabel("time")
    fig.suptitle("Lab 1")
    _save(fig, "lab1_introductory.png")


def plot_lab2():
    sol, _, p = lab2.solve(**lab2.CASE_A10)
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0], label="optimal")
    ax[0].plot(sol.time, lab2.untreated_state(sol.time, p["r"], p["M"], p["x0"]), "k--", label=r"$u\equiv 0$")
    ax[0].set_ylabel("mold $x$")
    ax[0].legend(fontsize=8)
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[1].set_ylabel("$u$")
    ax[1].set_xlabel("time")
    fig.suptitle("Lab 2 ($A=10$)")
    _save(fig, "lab2_mold_fungicide.png")


def plot_lab3():
    sol, _, _ = lab3.solve()
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0])
    ax[0].set_ylabel("bacteria $x$")
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[1].set_ylabel("$u$")
    ax[1].set_xlabel("time")
    fig.suptitle("Lab 3")
    _save(fig, "lab3_bacteria.png")


def plot_lab4():
    sol, _, _ = lab4.solve()
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0])
    ax[0].set_ylabel("$x$")
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[1].set_ylabel("$u$")
    ax[1].set_xlabel("time")
    fig.suptitle("Lab 4 (wide bounds)")
    _save(fig, "lab4_bounded.png")


def plot_lab5():
    sol, _, p = lab5.solve()
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0], label="treated")
    ax[0].plot(sol.time, lab5.untreated_N(sol.time, p["r"], p["N0"]), "k--", label=r"$u\equiv 0$")
    ax[0].set_ylabel("tumor $N$")
    ax[0].legend(fontsize=8)
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[1].set_ylabel("$u$")
    ax[1].set_xlabel("time")
    fig.suptitle("Lab 5")
    _save(fig, "lab5_cancer.png")


def plot_lab6():
    sol, _, p = lab6.solve()
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0])
    ax[0].set_ylabel("fish $x$")
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[1].set_ylabel("harvest $u$")
    ax[1].set_xlabel("time")
    fig.suptitle("Lab 6")
    _save(fig, "lab6_fish.png")


def plot_lab7():
    sol, _, _ = lab7.solve()
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0], label="$S$")
    ax[0].plot(sol.time, sol.x[:, 2], label="$I$")
    ax[0].legend(fontsize=8)
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[1].set_ylabel("vaccination $u$")
    ax[1].set_xlabel("time")
    fig.suptitle("Lab 7")
    _save(fig, "lab7_epidemic.png")


def plot_lab8():
    sol, _, _ = lab8.solve()
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0], label="$T$")
    ax[0].plot(sol.time, sol.x[:, 2], label="$V$")
    ax[0].legend(fontsize=8)
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[1].set_ylabel("$u$ (0 = max therapy)")
    ax[1].set_xlabel("time")
    fig.suptitle("Lab 8")
    _save(fig, "lab8_hiv.png")


def plot_lab9():
    sol, _, _ = lab9.solve()
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0], label="park")
    ax[0].plot(sol.time, sol.x[:, 1], label="forest")
    ax[0].plot(sol.time, sol.x[:, 2], label="outside")
    ax[0].legend(fontsize=8)
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post", label="$u_p$")
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 1], where="post", label="$u_f$")
    ax[1].legend(fontsize=8)
    ax[1].set_xlabel("time")
    fig.suptitle("Lab 9")
    _save(fig, "lab9_bears.png")


def plot_lab10():
    sol, _, p = lab10.solve()
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0])
    ax[0].axhline(p["l"], color="k", ls="--")
    ax[0].set_ylabel("glucose $x_1$")
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[1].set_ylabel("insulin $u$")
    ax[1].set_xlabel("time (d)")
    fig.suptitle("Lab 10 ($T=0.2$)")
    _save(fig, "lab10_glucose.png")


def plot_lab11():
    sol, _, p = lab11.solve()
    tsw = lab11.analytic_switch(p["k"], p["r"], p["T"])
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0])
    ax[0].set_ylabel("timber $x$")
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[1].axvline(tsw, color="k", ls="--")
    ax[1].set_ylabel("$u$")
    ax[1].set_xlabel("time")
    fig.suptitle("Lab 11")
    _save(fig, "lab11_timber.png")


def plot_lab12():
    sol, _, p = lab12.solve(**lab12.INJECT)
    z = lab12.recover_z(sol, p["K"], p["z0"])
    fig, ax = plt.subplots(3, 1, figsize=(6.0, 6.2), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0])
    ax[0].set_ylabel("bacteria $x$")
    ax[1].plot(sol.time, z)
    ax[1].set_ylabel("contaminant $z$")
    ax[2].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[2].set_ylabel("$u$")
    ax[2].set_xlabel("time")
    fig.suptitle("Lab 12 ($K=2.1$)")
    _save(fig, "lab12_bioreactor.png")


def plot_lab13():
    sol, _, _ = lab13.solve()
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
    ax[0].plot(sol.time, sol.x[:, 0], label="prey")
    ax[0].plot(sol.time, sol.x[:, 1], label="predator")
    ax[0].legend(fontsize=8)
    ax[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[1].set_ylabel("pesticide $u$")
    ax[1].set_xlabel("time")
    fig.suptitle("Lab 13")
    _save(fig, "lab13_predator_prey.png")


if __name__ == "__main__":
    for fn in (
        plot_lab1, plot_lab2, plot_lab3, plot_lab4, plot_lab5, plot_lab6,
        plot_lab7, plot_lab8, plot_lab9, plot_lab10, plot_lab11, plot_lab12,
        plot_lab13,
    ):
        fn()
