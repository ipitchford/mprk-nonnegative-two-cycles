#!/usr/bin/env python3
"""EXPLORATORY (non-rigorous, float64): where does two-state MPRK22(alpha) have two-cycles?

For a grid of (alpha, K) at several equilibrium fractions p, count period-two
orbits of the scalar map F_{K,p} in (p,1) (sign changes of F^2 - id on a
fine log grid) and test whether one is attracting (|multiplier| < 1).
Writes X1_phase_diagram.csv and X1_phase_diagram.png.  Not a proof: grid-based
root detection can miss roots, and float64 cannot resolve cycle points closer
than ~1e-15 to 0 or 1.
"""
import csv, os, tempfile
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mprk22-mpl"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pds_schemes import scalar_F, two_cycles

HERE = Path(__file__).resolve().parent
alphas = np.round(np.arange(0.500, 0.6601, 0.004), 4)
logK = np.round(np.arange(1.0, 10.01, 0.25), 3)
ps = [1e-2, 1e-3, 1e-4, 1e-6]
rows = []
for p in ps:
    for a in alphas:
        for lk in logK:
            cyc = two_cycles(scalar_F(a, 10.0 ** lk, p), p)
            att = sum(1 for c in cyc if abs(c[2]) < 1)
            rows.append((p, a, lk, len(cyc), att))
with open(HERE / "X1_phase_diagram.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["p", "alpha", "log10_K", "n_two_cycles_found", "n_attracting"])
    w.writerows(rows)

plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})
fig, axes = plt.subplots(1, 4, figsize=(12, 3.6), sharey=True, constrained_layout=True)
cmap = matplotlib.colors.ListedColormap(["#f2f2f2", "#9ecae1", "#08519c"])
for ax, p in zip(axes, ps):
    Z = np.zeros((len(logK), len(alphas)))
    for (pp, a, lk, n, att) in rows:
        if pp != p:
            continue
        i = int(np.argmin(abs(logK - lk))); j = int(np.argmin(abs(alphas - a)))
        Z[i, j] = 2 if att > 0 else (1 if n > 0 else 0)
    da, dk = alphas[1] - alphas[0], logK[1] - logK[0]
    ax.imshow(Z, origin="lower", aspect="auto", cmap=cmap, vmin=0, vmax=2,
              extent=[alphas[0] - da / 2, alphas[-1] + da / 2, logK[0] - dk / 2, logK[-1] + dk / 2])
    ax.axvline((np.sqrt(5) - 1) / 2, color="k", lw=0.8, ls="--")
    ax.axvline(4 / 7, color="#a50f15", lw=0.8, ls=":")
    ax.set_title("p = %g" % p)
    ax.set_xlabel(r"$\alpha$")
axes[0].set_ylabel(r"$\log_{10} K$")
handles = [matplotlib.patches.Patch(color="#f2f2f2", label="no two-cycle found"),
           matplotlib.patches.Patch(color="#9ecae1", label="two-cycles, none attracting"),
           matplotlib.patches.Patch(color="#08519c", label="attracting two-cycle"),
           matplotlib.lines.Line2D([], [], color="k", ls="--", lw=0.8, label=r"$\alpha=1/\varphi$"),
           matplotlib.lines.Line2D([], [], color="#a50f15", ls=":", lw=0.8, label=r"$\alpha=4/7$")]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=5, frameon=False)
fig.savefig(HERE / "X1_phase_diagram.png", dpi=200, bbox_inches="tight")
# symmetric equilibrium p = 1/2 (the setting of the original experiments)
sym = 0
for a in alphas:
    for lk in logK:
        sym += len(two_cycles(scalar_F(a, 10.0 ** lk, 0.5), 0.5))
(HERE / "X1_symmetric_p_half.txt").write_text(
    "two-cycles found at p = 1/2 over the same (alpha, log10 K) grid: %d\n" % sym)
print("wrote X1_phase_diagram.csv/png;", len(rows), "parameter points; p=1/2 two-cycles:", sym)
