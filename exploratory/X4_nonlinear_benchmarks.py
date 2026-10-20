#!/usr/bin/env python3
"""EXPLORATORY (non-rigorous, float64): fixed large steps on nonlinear PDS.

(a) Reversible dimerisation written in atom variables, x = (a, z), z = 2[B]:
    production of z from a at rate k1 a^2, of a from z at rate k2 z; a + z = 1.
    Equilibrium k1 a^2 = k2 z; k1/k2 = 1e4 gives a* ~ 1e-2, 1e6 gives a* ~ 1e-3.
(b) Robertson's stiff chemical kinetics (conservative; equilibrium (0,0,1)).
Fixed-step MPRK22(alpha) for 20000 steps from the stated initial state; the
last iterates are classified as converged / two-cycle / other.
Writes X4_nonlinear_benchmarks.csv.
"""
import csv
from pathlib import Path
import numpy as np
from pds_schemes import mprk22

HERE = Path(__file__).resolve().parent


def dimer(k1, k2):
    def P(y):
        a, z = y
        return np.array([[0.0, k2 * z], [k1 * a * a, 0.0]])
    return P


def robertson(y):
    y1, y2, y3 = y
    P = np.zeros((3, 3))
    P[1, 0] = 0.04 * y1
    P[0, 1] = 1e4 * y2 * y3
    P[2, 1] = 3e7 * y2 * y2
    return P


def classify(Pf, y0, dt, alpha, steps=20000):
    y = y0.copy()
    hist = []
    for n in range(steps):
        y = mprk22(Pf, y, dt, alpha)
        if n >= steps - 3:
            hist.append(y.copy())
    y0_, y1_, y2_ = hist
    d1 = np.max(np.abs(y2_ - y1_) / np.maximum(np.abs(y2_), 1e-300))
    d2 = np.max(np.abs(y2_ - y0_) / np.maximum(np.abs(y2_), 1e-300))
    if d1 < 1e-10:
        return "fixed", y2_
    if d2 < 1e-10:
        return "two-cycle", y2_
    return "other(d1=%.1e,d2=%.1e)" % (d1, d2), y2_


def main():
    rows = []
    cases = [("dimer_k1/k2=1e4", dimer(1e4, 1.0), np.array([0.5, 0.5])),
             ("dimer_k1/k2=1e6", dimer(1e6, 1.0), np.array([0.5, 0.5])),
             ("robertson", robertson, np.array([1.0, 1e-20, 1e-20]))]
    for name, Pf, y0 in cases:
        for alpha in [0.52, 0.55, 4 / 7, 0.6, 2 / 3, 1.0]:
            for dt in [1.0, 1e2, 1e4, 1e6]:
                c, y = classify(Pf, y0, dt, alpha)
                rows.append((name, round(alpha, 6), dt, c, " ".join("%.4e" % v for v in y)))
                print(name, "alpha=%.4f dt=%g -> %s  %s" % (alpha, dt, c, rows[-1][-1]), flush=True)
    with open(HERE / "X4_nonlinear_benchmarks.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "alpha", "dt", "outcome", "final_state"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
