#!/usr/bin/env python3
"""EXPLORATORY (non-rigorous; mpmath and float64): supporting numerical checks
quoted in the paper.  Writes X6_supporting_checks.json.

(a) basin sample for the three-species chain of Theorem 9 (200 random starts);
(b) one-step ratios of the log-odds distance to equilibrium, two species, by alpha;
(c) one-step ratios of Hilbert's projective distance to equilibrium, three species;
(d) sign of the Schwarzian derivative of the witness map F on a log grid;
(e) local multiplier at the equilibrium of the dimerisation test (script X4);
(f) observed convergence orders of the implemented schemes.
"""
import json, random
from pathlib import Path
import numpy as np
from mpmath import mp
from pds_schemes import mprk22, mprk43, sspmprk22

HERE = Path(__file__).resolve().parent
out = {}
mp.dps = 50


def mprk22_mp(y, A, h, a):
    n = len(y); r = 1 / a; b = 1 - 1 / (2 * a); g = 1 / (2 * a); I = mp.eye(n)
    u = mp.lu_solve(I - a * h * A, y)
    q = [(b * y[i] + g * u[i]) / (y[i] ** (1 - r) * u[i] ** r) for i in range(n)]
    return mp.lu_solve(I - h * A * mp.diag(q), y)


# (a)
A = mp.matrix([[-10, 10000, 0], [10, -20000, 1], [0, 10000, -1]])
a47 = mp.mpf(4) / 7
pi = mp.matrix([mp.mpf(1000) / 11001, mp.mpf(1) / 11001, mp.mpf(10000) / 11001])
random.seed(1)
cyc = eq = other = 0
mp.dps = 30
for s in range(200):
    w = [random.random() for _ in range(3)]
    y = mp.matrix([v / sum(w) for v in w])
    for n in range(600):
        y = mprk22_mp(y, A, 1, a47)
    y2 = mprk22_mp(mprk22_mp(y, A, 1, a47), A, 1, a47)
    if mp.norm(y - pi) < 1e-15:
        eq += 1
    elif mp.norm(y2 - y) < 1e-15:
        cyc += 1
    else:
        other += 1
out["a_three_species_basin_sample"] = {"to_equilibrium": eq, "to_two_cycle": cyc, "other": other}
print("(a)", out["a_three_species_basin_sample"], flush=True)

# (b)
mp.dps = 50
def F2s(t, a, K, p):
    r = 1 / a; b = 1 - 1 / (2 * a); g = 1 / (2 * a); u = (t + a * K * p) / (1 + a * K)
    q1 = (b * t + g * u) / (t ** (1 - r) * u ** r)
    q2 = (b * (1 - t) + g * (1 - u)) / ((1 - t) ** (1 - r) * (1 - u) ** r)
    return (t + K * p * q2) / (1 + K * ((1 - p) * q1 + p * q2))
lo = lambda x: mp.log(x / (1 - x))
random.seed(3)
res = {}
for a in ["1/2", "0.55", "4/7", "0.6", "0.62", "2/3", "0.75", "0.9", "1", "1.25", "1.5", "2", "3"]:
    A_ = mp.mpf(eval(a)) if "/" in a else mp.mpf(a)
    worst = 0
    for _ in range(3000):
        p = mp.mpf(10) ** random.uniform(-8, -0.3); K = mp.mpf(10) ** random.uniform(-1, 9)
        t = mp.mpf(10) ** random.uniform(-25, -0.3) if random.random() < 0.5 else 1 - mp.mpf(10) ** random.uniform(-25, -0.3)
        if abs(t - p) < 1e-40:
            continue
        worst = max(worst, abs(lo(F2s(t, A_, K, p)) - lo(p)) / abs(lo(t) - lo(p)))
    res[a] = float(worst)
out["b_max_logodds_ratio_two_species"] = res
print("(b)", res, flush=True)

# (c)
def stat(A):
    n = A.rows; M = A.copy()
    for j in range(n):
        M[n - 1, j] = 1
    return mp.lu_solve(M, mp.matrix([0] * (n - 1) + [1]))
def dH(y, p):
    v = [y[i] / p[i] for i in range(len(p))]
    return mp.log(max(v) / min(v))
random.seed(11)
res = {}
for a in [mp.mpf(1) / 2, mp.mpf(1)]:
    worst = 0
    for trial in range(150):
        A_ = mp.matrix(3, 3)
        for i in range(3):
            for j in range(3):
                if i != j:
                    A_[i, j] = mp.mpf(10) ** random.uniform(-5, 5)
        for j in range(3):
            A_[j, j] = -sum(A_[i, j] for i in range(3) if i != j)
        p_ = stat(A_); h = mp.mpf(10) ** random.uniform(-3, 8)
        for s in range(8):
            y = mp.matrix([mp.mpf(10) ** random.uniform(-15, 0) for i in range(3)]); y = y / sum(y)
            d0 = dH(y, p_)
            if d0 < 1e-20:
                continue
            worst = max(worst, dH(mprk22_mp(y, A_, h, a), p_) / d0)
    res[str(float(a))] = float(worst)
out["c_max_hilbert_ratio_three_species"] = res
print("(c)", res, flush=True)

# (d)
Fw = lambda t: F2s(t, mp.mpf(4) / 7, 10000, mp.mpf(1) / 1000)
def S(t):
    d1 = mp.diff(Fw, t, 1); d2 = mp.diff(Fw, t, 2); d3 = mp.diff(Fw, t, 3)
    return d3 / d1 - mp.mpf(3) / 2 * (d2 / d1) ** 2
pts = [mp.mpf(10) ** (-12 + 12 * i / 200) for i in range(201)]
pts = pts + [1 - x for x in pts if x < 0.5]
vals = [S(x) for x in pts if 0 < x < 1]
out["d_schwarzian_witness"] = {"grid_points": len(vals), "positive": sum(1 for v in vals if v > 0),
                               "max": float(max(vals)), "min": float(min(vals))}
print("(d)", out["d_schwarzian_witness"], flush=True)

# (e)
def dimer_step(y, dt, al, k1, k2):
    def P(y):
        a, z = y; return [[0, k2 * z], [k1 * a * a, 0]]
    def pat(v, Pm, sig):
        M = mp.matrix([[1 + dt * Pm[1][0] / sig[0], -dt * Pm[0][1] / sig[1]],
                       [-dt * Pm[1][0] / sig[0], 1 + dt * Pm[0][1] / sig[1]]])
        return mp.lu_solve(M, mp.matrix(v))
    P0 = P(y); u = pat(y, [[al * x for x in r] for r in P0], y)
    r = 1 / al; sig = [y[i] ** (1 - r) * u[i] ** r for i in range(2)]
    P1 = P(u); b1 = 1 - 1 / (2 * al); b2 = 1 / (2 * al)
    return pat(y, [[b1 * P0[i][j] + b2 * P1[i][j] for j in range(2)] for i in range(2)], sig)
k1, k2 = mp.mpf(10) ** 4, mp.mpf(1)
astar = (-k2 + mp.sqrt(k2 ** 2 + 4 * k1 * k2)) / (2 * k1)
res = {}
for al in ["0.52", "0.55", "4/7", "0.6", "2/3", "1"]:
    A_ = mp.mpf(eval(al)) if "/" in al else mp.mpf(al)
    for dt in [1, 100]:
        f = lambda x: dimer_step(mp.matrix([x, 1 - x]), mp.mpf(dt), A_, k1, k2)[0]
        res["alpha=%s,dt=%d" % (al, dt)] = float(mp.diff(f, astar))
out["e_dimer_local_multiplier"] = res
print("(e)", res, flush=True)

# (f)
def Pd(y):
    a, z = y; return np.array([[0.0, 1.0 * z], [2.0 * a * a, 0.0]])
def ref():
    y = np.array([0.9, 0.1]); n = 20000
    for i in range(n):
        y = mprk43(Pd, y, 1.0 / n, 0.74)
    return y
yr = ref()
orders = {}
for name, st in [("MPRK43II(0.74)", lambda y, h: mprk43(Pd, y, h, 0.74)),
                 ("MPRK22(4/7)", lambda y, h: mprk22(Pd, y, h, 4 / 7)),
                 ("SSPMPRK22(0.5,1)", lambda y, h: sspmprk22(Pd, y, h, 0.5, 1.0))]:
    errs = []
    for n in [10, 20, 40, 80]:
        y = np.array([0.9, 0.1])
        for i in range(n):
            y = st(y, 1.0 / n)
        errs.append(float(np.abs(y - yr).max()))
    orders[name] = [round(float(np.log2(errs[i] / errs[i + 1])), 2) for i in range(3)]
out["f_observed_orders"] = orders
print("(f)", orders, flush=True)
(HERE / "X6_supporting_checks.json").write_text(json.dumps(out, indent=2) + "\n")
