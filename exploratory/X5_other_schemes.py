#!/usr/bin/env python3
"""EXPLORATORY (non-rigorous, float64): do other positive-tableau Patankar schemes
show the same asymmetric large-step two-cycles?

Two-state linear test y' = K A_p y (dt = 1).  For each scheme the scalar
first-component map is built from the full scheme (vectorised 2x2 Patankar
solves) and period-two points are located from sign changes of F^2 - id.
Schemes (formulas as in pds_schemes.py): MPRK22(alpha) (reference),
SSPMPRK22(a21, b10) with a21 > 0, and MPRK43II(gamma).
Writes X5_other_schemes.csv.
"""
import csv
from pathlib import Path
import numpy as np
from pds_schemes import two_cycles

HERE = Path(__file__).resolve().parent


def pat2(v1, v2, P21, P12, s1, s2):
    a, b = P21 / s1, P12 / s2          # dt = 1
    det = 1 + a + b
    return ((1 + b) * v1 + b * v2) / det, (a * v1 + (1 + a) * v2) / det


def prod(K, p, y1, y2):
    return K * (1 - p) * y1, K * p * y2   # P21 (2 from 1), P12 (1 from 2)


def F_mprk22(alpha, K, p):
    def F(t):
        y1, y2 = t, 1 - t
        A21, A12 = prod(K, p, y1, y2)
        u1, u2 = pat2(y1, y2, alpha * A21, alpha * A12, y1, y2)
        r = 1 / alpha
        s1, s2 = y1 ** (1 - r) * u1 ** r, y2 ** (1 - r) * u2 ** r
        B21, B12 = prod(K, p, u1, u2)
        b1, b2 = 1 - 1 / (2 * alpha), 1 / (2 * alpha)
        return pat2(y1, y2, b1 * A21 + b2 * B21, b1 * A12 + b2 * B12, s1, s2)[0]
    return F


def F_ssp(a21, b10, K, p):
    b20 = 1 - 1 / (2 * b10) - a21 * b10
    b21 = 1 / (2 * b10)
    s = (b20 + b21 + a21 * b10 ** 2) / (b10 * (b20 + b21))
    def F(t):
        y1, y2 = t, 1 - t
        A21, A12 = prod(K, p, y1, y2)
        u1, u2 = pat2(y1, y2, b10 * A21, b10 * A12, y1, y2)
        s1, s2 = y1 ** (1 - s) * u1 ** s, y2 ** (1 - s) * u2 ** s
        B21, B12 = prod(K, p, u1, u2)
        v1, v2 = (1 - a21) * y1 + a21 * u1, (1 - a21) * y2 + a21 * u2
        return pat2(v1, v2, b20 * A21 + b21 * B21, b20 * A12 + b21 * B12, s1, s2)[0]
    return F


def F_mprk43(gamma, K, p):
    a21 = 2.0 / 3; a31 = a21 - 1 / (4 * gamma); a32 = 1 / (4 * gamma)
    b1, b2, b3 = 0.25, 0.75 - gamma, gamma
    beta2 = 1 / (2 * a21); beta1 = 1 - beta2
    q1 = 1 / (3 * a21 * (a31 + a32) * b3); q2 = 1 / a21
    def F(t):
        y1, y2 = t, 1 - t
        A21, A12 = prod(K, p, y1, y2)
        u1, u2 = pat2(y1, y2, a21 * A21, a21 * A12, y1, y2)
        s1, s2 = y1 ** (1 - q1) * u1 ** q1, y2 ** (1 - q1) * u2 ** q1
        B21, B12 = prod(K, p, u1, u2)
        w1, w2 = pat2(y1, y2, a31 * A21 + a32 * B21, a31 * A12 + a32 * B12, s1, s2)
        r1, r2 = y1 ** (1 - q2) * u1 ** q2, y2 ** (1 - q2) * u2 ** q2
        g1, g2 = pat2(y1, y2, beta1 * A21 + beta2 * B21, beta1 * A12 + beta2 * B12, r1, r2)
        C21, C12 = prod(K, p, w1, w2)
        return pat2(y1, y2, b1 * A21 + b2 * B21 + b3 * C21, b1 * A12 + b2 * B12 + b3 * C12, g1, g2)[0]
    return F


def main():
    schemes = [("MPRK22(4/7)", lambda K, p: F_mprk22(4 / 7, K, p)),
               ("MPRK22(1)", lambda K, p: F_mprk22(1.0, K, p)),
               ("SSPMPRK22(0.5,1)", lambda K, p: F_ssp(0.5, 1.0, K, p)),
               ("SSPMPRK22(0.2,0.6)", lambda K, p: F_ssp(0.2, 0.6, K, p)),
               ("SSPMPRK22(0.1,0.55)", lambda K, p: F_ssp(0.1, 0.55, K, p)),
               ("MPRK43II(0.5)", lambda K, p: F_mprk43(0.5, K, p)),
               ("MPRK43II(0.563)", lambda K, p: F_mprk43(0.563, K, p)),
               ("MPRK43II(0.75)", lambda K, p: F_mprk43(0.75, K, p)),
               ("MPRK43II(0.375)", lambda K, p: F_mprk43(0.375, K, p))]
    rows = []
    for name, mk in schemes:
        for p in [0.5, 1e-2, 1e-3, 1e-6]:
            res = []
            for lk in [1, 2, 3, 4, 5, 6, 8]:
                with np.errstate(all="ignore"):
                    c = two_cycles(mk(10.0 ** lk, p), p)
                att = sum(1 for x in c if abs(x[2]) < 1)
                res.append("%d/%d" % (att, len(c)))
                rows.append((name, p, lk, len(c), att))
            print("%-20s p=%-6g attracting/all two-cycles at log10K=1,2,3,4,5,6,8: %s" % (name, p, res), flush=True)
    with open(HERE / "X5_other_schemes.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["scheme", "p", "log10_K", "n_two_cycles_found", "n_attracting"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
