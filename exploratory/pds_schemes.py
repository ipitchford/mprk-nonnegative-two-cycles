"""EXPLORATORY NUMERICS ONLY (floating point; not part of any proof).

Modified Patankar schemes for conservative production-destruction systems
y_i' = sum_j (P_ij(y) - P_ji(y)), with P_ij >= 0 the rate of production of
species i from species j.  Implemented following the formulas of the cited
papers, cross-checked against the PositiveIntegrators.jl implementation
(src/mprk.jl, src/sspmprk.jl):

  MPRK22(alpha)        Kopecz & Meister (2018), Appl. Numer. Math. 123
  SSPMPRK22(a, b)      Huang & Shu (2019), J. Sci. Comput. 78
  MPRK43II(gamma)      Kopecz & Meister (2018), BIT 58

Embedded error estimate of MPRK22 (Kopecz & Meister 2018; PositiveIntegrators.jl):
y_emb = sigma = y^(1-1/alpha) * u^(1/alpha), error = y+ - sigma.
Works with numpy arrays (float64) or mpmath matrices via the `backend`.
"""
import numpy as np


def patankar_step(v, P, sigma, dt):
    """Solve M y = v, M_ii = 1 + dt sum_j P_ji / sigma_i, M_ij = -dt P_ij / sigma_j."""
    n = len(v)
    Mtx = -dt * P / sigma[None, :]
    np.fill_diagonal(Mtx, 0.0)
    destr = P.sum(axis=0) - np.diag(P)          # sum_j P_ji over j != i
    Mtx[np.diag_indices(n)] = 1.0 + dt * destr / sigma
    return np.linalg.solve(Mtx, v)


def linear_P(A):
    """Production matrix function for y' = A y (A conservative Metzler)."""
    Aoff = A.copy()
    np.fill_diagonal(Aoff, 0.0)
    return lambda y: Aoff * y[None, :]


def mprk22(Pf, y, dt, alpha, embedded=False):
    P0 = Pf(y)
    u = patankar_step(y, alpha * P0, y, dt)
    r = 1.0 / alpha
    sigma = y ** (1 - r) * u ** r
    P1 = Pf(u)
    yp = patankar_step(y, (1 - 1 / (2 * alpha)) * P0 + 1 / (2 * alpha) * P1, sigma, dt)
    if embedded:
        return yp, sigma
    return yp


def sspmprk22(Pf, y, dt, a21, b10):
    a10, a20 = 1.0, 1 - a21
    b20 = 1 - 1 / (2 * b10) - a21 * b10
    b21 = 1 / (2 * b10)
    s = (b20 + b21 + a21 * b10 ** 2) / (b10 * (b20 + b21))
    P0 = Pf(y)
    u = patankar_step(a10 * y, b10 * P0, y, dt)
    sigma = y ** (1 - s) * u ** s
    P1 = Pf(u)
    return patankar_step(a20 * y + a21 * u, b20 * P0 + b21 * P1, sigma, dt)


def mprk43(Pf, y, dt, gamma):
    a21 = 2.0 / 3
    a31 = a21 - 1 / (4 * gamma)
    a32 = 1 / (4 * gamma)
    b1, b2, b3 = 0.25, 0.75 - gamma, gamma
    beta2 = 1 / (2 * a21)
    beta1 = 1 - beta2
    q1 = 1 / (3 * a21 * (a31 + a32) * b3)
    q2 = 1 / a21
    P0 = Pf(y)
    u2 = patankar_step(y, a21 * P0, y, dt)
    sig = y ** (1 - q1) * u2 ** q1
    P1 = Pf(u2)
    u3 = patankar_step(y, a31 * P0 + a32 * P1, sig, dt)
    sig2 = y ** (1 - q2) * u2 ** q2
    sig3 = patankar_step(y, beta1 * P0 + beta2 * P1, sig2, dt)
    P2 = Pf(u3)
    return patankar_step(y, b1 * P0 + b2 * P1 + b3 * P2, sig3, dt)


def scalar_F(alpha, K, p):
    """Two-state MPRK22 first-component map (paper eq. for F_{K,p}), vectorised in t."""
    r = 1.0 / alpha
    be, ga = 1 - 1 / (2 * alpha), 1 / (2 * alpha)
    kap = 0.0 if np.isinf(K) else 1.0 / K

    def F(t):
        u = (kap * t + alpha * p) / (kap + alpha)
        omu = (kap * (1 - t) + alpha * (1 - p)) / (kap + alpha)
        v1, v2 = t / u, (1 - t) / omu
        q1 = v1 ** (r - 1) * (be * v1 + ga)
        q2 = v2 ** (r - 1) * (be * v2 + ga)
        return (kap * t + p * q2) / (kap + (1 - p) * q1 + p * q2)
    return F


def two_cycles(F, p, n=6000):
    """Approximate period-two points of F in (p,1) via sign changes of F(F(t)) - t.

    Returns a list of (t, F(t), multiplier).  Floating point; grid-based."""
    x = np.concatenate([p + p * np.logspace(-10, np.log10((0.5 - p) / p), n // 2),
                        1 - np.logspace(-15, np.log10(0.5), n // 2)[::-1]])
    x = np.unique(x[(x > p) & (x < 1)])
    with np.errstate(all="ignore"):
        G = F(F(x)) - x
    s = np.sign(G)
    idx = np.where((s[:-1] * s[1:]) < 0)[0]
    out = []
    for i in idx:
        lo, hi = x[i], x[i + 1]
        slo = s[i]
        for _ in range(80):
            m = 0.5 * (lo + hi)
            with np.errstate(all="ignore"):
                gm = F(F(m)) - m
            if np.sign(gm) == slo:
                lo = m
            else:
                hi = m
        t = 0.5 * (lo + hi)
        ft = F(t)
        # discard rounding-noise sign changes at the fixed point itself
        if not (t > p * (1 + 1e-4) and ft < p * (1 - 1e-4)):
            continue
        def d(y):
            hh = min(y, 1 - y) * 1e-6
            return (F(y + hh) - F(y - hh)) / (2 * hh)
        out.append((t, ft, d(t) * d(ft)))
    return out
