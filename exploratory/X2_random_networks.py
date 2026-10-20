#!/usr/bin/env python3
"""EXPLORATORY (non-rigorous, float64): how often does MPRK22(alpha) end on a
two-cycle for random stiff linear networks?

N-species conservative linear systems y' = A y with every off-diagonal rate
log-uniform in [1e-3, 1e3] (fixed seed), fixed step h, four initial states
(uniform and three random points of the simplex), 5000 steps.  Final state
classified as: equilibrium (relative error < 1e-9), two-cycle (y_{n+2} = y_n to
1e-11 relative but y_{n+1} != y_n), or unsettled.  Iterates are renormalised to unit mass after every step (the exact scheme
conserves mass; without this, rounding in the column sums of A causes a
spurious drift when h is large).  Unsettled runs are further split by whether the error to the equilibrium was
at least halved over the last 1000 steps (slow convergence) or not.
Writes X2_random_networks.csv.
"""
import csv
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
rng = np.random.default_rng(20261020)


def make_systems(N, count):
    As = []
    for _ in range(count):
        R = 10.0 ** rng.uniform(-3, 3, size=(N, N))
        np.fill_diagonal(R, 0.0)
        A = R - np.diag(R.sum(axis=0))
        As.append(A)
    return np.array(As)


def stationary(A):
    N = A.shape[0]
    M = A.copy(); M[-1, :] = 1.0
    b = np.zeros(N); b[-1] = 1.0
    return np.linalg.solve(M, b)


def run(As, Y0, h, alpha, steps=5000):
    B, N, _ = As.shape
    I = np.eye(N)[None]
    S = I - alpha * h * As
    r, be, ga = 1 / alpha, 1 - 1 / (2 * alpha), 1 / (2 * alpha)
    y = Y0.copy()
    hist = []
    for n in range(steps):
        u = np.linalg.solve(S, y[..., None])[..., 0]
        v = y / u
        q = v ** (r - 1) * (be * v + ga)
        Mt = I - h * As * q[:, None, :]
        y = np.linalg.solve(Mt, y[..., None])[..., 0]
        y = np.maximum(y, 1e-300)
        # MPRK22 conserves mass exactly; renormalise to remove floating-point drift
        # caused by rounding in the column sums of A (amplified by large h)
        y = y / y.sum(axis=1, keepdims=True)
        if n == steps - 1001:
            early = y.copy()
        if n >= steps - 3:
            hist.append(y.copy())
    return hist + [early]  # y_{n-2}, y_{n-1}, y_n, y_{n-1000}


def main():
    rows = []
    summary = []
    for N, count in [(3, 100), (4, 50)]:
        As = make_systems(N, count)
        pis = np.array([stationary(A) for A in As])
        starts = [np.full(N, 1.0 / N)] + [rng.dirichlet(np.ones(N)) for _ in range(3)]
        for h in [10.0, 1e3, 1e5]:
            for alpha in [0.5, 0.52, 0.55, 4 / 7, 0.6, 0.62, 2 / 3, 0.8, 1.0]:
                cls_sys = np.zeros((count, 4), dtype=int)
                for y0 in starts:
                    Y0 = np.repeat(y0[None], count, axis=0)
                    y0_, y1_, y2_, ye_ = run(As, Y0, h, alpha)
                    err_e = np.linalg.norm(ye_ - pis, axis=1) / np.linalg.norm(pis, axis=1)
                    err = np.linalg.norm(y2_ - pis, axis=1) / np.linalg.norm(pis, axis=1)
                    per2 = np.linalg.norm(y2_ - y0_, axis=1) / np.linalg.norm(y2_, axis=1)
                    per1 = np.linalg.norm(y2_ - y1_, axis=1) / np.linalg.norm(y2_, axis=1)
                    eq = err < 1e-9
                    cyc = (~eq) & (per2 < 1e-11) & (per1 > 1e-6)
                    uns = ~(eq | cyc)
                    slow = uns & (err < 0.5 * err_e)
                    cls_sys[:, 0] += eq; cls_sys[:, 1] += cyc; cls_sys[:, 2] += uns; cls_sys[:, 3] += slow
                any_cycle = int((cls_sys[:, 1] > 0).sum())
                all_eq = int((cls_sys[:, 0] == len(starts)).sum())
                any_uns = int((cls_sys[:, 2] > 0).sum())
                uns_not_decr = int(((cls_sys[:, 2] - cls_sys[:, 3]) > 0).sum())
                rows.append((N, count, h, round(alpha, 6), all_eq, any_cycle, any_uns, uns_not_decr))
                print("N=%d h=%g alpha=%.4f: all starts->eq %d/%d, some start->2-cycle %d, some unsettled %d (of which error not halved in last 1000 steps: %d)"
                      % (N, h, alpha, all_eq, count, any_cycle, any_uns, uns_not_decr), flush=True)
    with open(HERE / "X2_random_networks.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["N", "systems", "h", "alpha", "systems_all_starts_to_equilibrium",
                    "systems_with_a_start_ending_on_two_cycle", "systems_with_unsettled_start",
                    "systems_with_unsettled_start_error_not_halved_in_last_1000_steps"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
