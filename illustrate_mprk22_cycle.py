#!/usr/bin/env python3
"""Illustrations only (NumPy/Matplotlib, decimal arithmetic); the proofs are the
exact verifiers.  Produces trajectory.csv, fig_trajectory.png, fig_portrait.png."""
from pathlib import Path
from decimal import Decimal as D, getcontext
import csv, os, tempfile
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mprk22-matplotlib-cache"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

getcontext().prec = 80
out = Path(__file__).parent


def step(t):
    u = (7 * t + 40) / 40007
    def q(r):
        return (r + 7) * (r * r * r).sqrt().sqrt() / 8
    q1, q2 = q(t / u), q((1 - t) / (1 - u))
    return (t + 10 * q2) / (1 + 9990 * q1 + 10 * q2)


t = D(1) / 2
rows = []
for n in range(41):
    rows.append((n, str(t)))
    t = step(t)
with (out / "trajectory.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["step", "first_component_80_decimal_precision"])
    w.writerows(rows)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                     "axes.spines.top": False, "axes.spines.right": False})
ns = np.array([n for n, _ in rows]); ys = np.array([float(v) for _, v in rows])
fig, ax = plt.subplots(figsize=(8, 4.2), constrained_layout=True)
ax.axhspan(0.48, 0.507, color="#76466F", alpha=0.15)
ax.axhspan(1e-7, 1.17e-7, color="#A5714B", alpha=0.25)
ax.text(40.6, 0.49, "H", va="center", fontsize=11)
ax.text(40.6, 1.08e-7, "L", va="center", fontsize=11)
ax.plot(ns, ys, color="#76466F", lw=1.2, marker="o", ms=3.2, label="MPRK22(4/7), step size h = 1")
tt = np.linspace(0, 40, 400)
ax.plot(tt, 1e-3 + 0.499 * np.exp(-1e4 * tt), color="#4C7A2F", lw=1.8, label="exact solution")
ax.axhline(1e-3, color="#4C7A2F", lw=0.8, ls=":")
ax.text(1.0, 2.6e-3, "equilibrium p = 0.001", color="#4C7A2F", fontsize=9,
        bbox=dict(facecolor="white", edgecolor="none", pad=1))
ax.set_yscale("log"); ax.set_ylim(5e-8, 1.2); ax.set_xlim(-0.5, 42)
ax.set_xlabel("time step n"); ax.set_ylabel("first component $y_1^n$")
ax.set_xticks([0, 10, 20, 30, 40]); ax.set_yticks([1e-7, 1e-5, 1e-3, 1e-1])
ax.grid(axis="y", alpha=0.18)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=10)
fig.savefig(out / "fig_trajectory.png", dpi=220)

# phase portrait of F^2 (floating point is adequate at these magnitudes)
def F(t):
    u = (7 * t + 40) / 40007
    W = lambda v: (v + 7) / 8 * v ** 0.75
    return (t + 10 * W((1 - t) / (1 - u))) / (1 + 9990 * W(t / u) + 10 * W((1 - t) / (1 - u)))
x = np.logspace(-8.5, np.log10(0.95), 3000)
G = F(F(x))
roots = {"t_L": 1.0814760e-7, "s_L": 1.9651694e-6, "p": 1e-3, "s_H": 0.10805097, "t_H": 0.49352223}
fig, ax = plt.subplots(figsize=(8, 4.6), constrained_layout=True)
ax.axvspan(roots["s_L"], roots["s_H"], color="#4C7A2F", alpha=0.12)
ax.plot(x, G / x, color="#76466F", lw=1.5, label=r"$F^2(t)/t$")
ax.axhline(1, color="k", lw=0.8)
for name, r in roots.items():
    ax.axvline(r, color="0.5", lw=0.6, ls="--")
    lab = {"t_L": r"$t_L$", "s_L": r"$s_L$", "p": r"$p$", "s_H": r"$s_H$", "t_H": r"$t_H$"}[name]
    ax.text(r * 1.1, 8.0, lab, fontsize=10)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("first component t"); ax.set_ylabel(r"$F^2(t)/t$")
ax.set_ylim(0.3, 12)
ax.text(2.2e-5, 0.42, "basin of the\nequilibrium\n" + r"$(s_L, s_H)$", color="#2f5a1c", fontsize=9, ha="center")
ax.text(1.0e-8 * 3, 0.42, "to the\ntwo-cycle", fontsize=9, ha="center")
ax.text(0.27, 0.42, "to the\ntwo-cycle", fontsize=9, ha="center")
ax.legend(loc="upper left", frameon=False)
fig.savefig(out / "fig_portrait.png", dpi=220)
print("Created trajectory.csv, fig_trajectory.png, fig_portrait.png (illustrations only).")
