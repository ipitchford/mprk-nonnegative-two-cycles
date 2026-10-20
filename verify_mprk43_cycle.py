#!/usr/bin/env python3
"""Certificate C7: an attracting two-cycle of the third-order scheme MPRK43II(37/50).

Scheme: MPRK43II(gamma) of Kopecz & Meister (BIT 58, 2018; arXiv:1703.05052,
Lemma 2.4 case II, Theorem 2.7 / eqs. (63), scheme (68) with delta = 1), as
written out in the paper (Section 6.3), with gamma = 37/50.  Coefficients
    a21 = 2/3, a31 = 2/3 - 1/(4 gamma) = 73/222, a32 = 1/(4 gamma) = 25/74,
    b1 = 1/4, b2 = 3/4 - gamma = 1/100, b3 = gamma = 37/50,
all nonzero coefficients positive.  Patankar denominators:
    rho  = y^(1-Q1) * y2^Q1   with Q1 = 1/(3 a21 (a31+a32) b3) = 75/74,
    mu   = y^(1-Q2) * y2^Q2   with Q2 = 1/a21 = 3/2,
    sigma = embedded MPRK22(a21) step with denominators mu, weights
            beta1 = 1 - 1/(2 a21), beta2 = 1/(2 a21),
and the final solve uses denominators sigma (function step() below).
The code also agrees with PositiveIntegrators.jl, src/mprk.jl at commit
37efc009f8baed988b999f56ee0dc76021bf8ee1 (cross-check only, not a dependency).

Test problem: y' = A y, A = [[-9999000, 1000], [9999000, -1000]], h = 1,
i.e. K = 10^7 and equilibrium fraction p = 10^-4 in the paper's notation.

The program evaluates the exact two-state step (all four linear 2x2 Patankar
solves written out) in outward-rounded rational interval arithmetic, together
with its t-derivative by forward-mode automatic differentiation on intervals,
and certifies for rational intervals L (near 0) and H (containing 1/2):
    F' < 0 on L and on H;  F(L) inside int H;  F(H) inside int L;
    sup_L |F'| * sup_H |F'| < 1,
and that the physical equilibrium p is a fixed point with |F'(p)| < 1.
Hence (as for the counterexample theorem) there is an attracting prime-period-two orbit
attracting all of L and H, and the initial state (1/2, 1/2) converges to it,
while the equilibrium is locally asymptotically stable.
Run:  python -S -O verify_mprk43_cycle.py --out NEW.json
"""
from fractions import Fraction as Fr
import mprk_exact as E
from mprk_exact import Iv, require, as_iv
import certlib

E.PREC = 160
GAMMA = Fr(37, 50)
A21 = Fr(2, 3)
A31 = A21 - 1 / (4 * GAMMA)
A32 = 1 / (4 * GAMMA)
B1, B2, B3 = Fr(1, 4), Fr(3, 4) - GAMMA, GAMMA
BETA2 = 1 / (2 * A21)
BETA1 = 1 - BETA2
Q1 = 1 / (3 * A21 * (A31 + A32) * B3)
Q2 = 1 / A21
K = Fr(10 ** 7)
P = Fr(1, 10 ** 4)


class D:
    """Forward-mode dual number over intervals: value and derivative."""
    __slots__ = ("v", "d")

    def __init__(self, v, d=0):
        self.v, self.d = as_iv(v), as_iv(d)

    @staticmethod
    def lift(x):
        return x if isinstance(x, D) else D(x, 0)

    def __add__(self, o):
        o = D.lift(o); return D(self.v + o.v, self.d + o.d)
    __radd__ = __add__

    def __neg__(self):
        return D(-self.v, -self.d)

    def __sub__(self, o):
        return self + (-D.lift(o))

    def __rsub__(self, o):
        return D.lift(o) + (-self)

    def __mul__(self, o):
        o = D.lift(o); return D(self.v * o.v, self.d * o.v + self.v * o.d)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = D.lift(o)
        q = self.v / o.v
        return D(q, (self.d - q * o.d) / o.v)

    def __rtruediv__(self, o):
        return D.lift(o) / self

    def rpow(self, e):
        """self ** e for a rational exponent e and a positive value."""
        e = Fr(e)
        k, m = e.numerator, e.denominator
        val = self.v.rpow(k, m)
        if e == 0:
            return D(1, 0)
        return D(val, e * val / self.v * self.d)


def rp(x, e):
    """x ** e for rational e on Iv or D."""
    if isinstance(x, D):
        return x.rpow(e)
    e = Fr(e)
    return as_iv(x).rpow(e.numerator, e.denominator)


def pat2(v1, v2, P21, P12, s1, s2):
    """Two-state Patankar solve with dt = 1."""
    a, b = P21 / s1, P12 / s2
    det = 1 + a + b
    return ((1 + b) * v1 + b * v2) / det, (a * v1 + (1 + a) * v2) / det


def prod(y1, y2):
    return K * (1 - P) * y1, K * P * y2      # P21 (2 from 1), P12 (1 from 2)


def step(t):
    """First component of one MPRK43II(gamma) step from (t, 1-t)."""
    y1, y2 = t, 1 - t
    A21_, A12_ = prod(y1, y2)
    u1, u2 = pat2(y1, y2, A21 * A21_, A21 * A12_, y1, y2)
    s1 = rp(y1, 1 - Q1) * rp(u1, Q1)
    s2 = rp(y2, 1 - Q1) * rp(u2, Q1)
    B21_, B12_ = prod(u1, u2)
    w1, w2 = pat2(y1, y2, A31 * A21_ + A32 * B21_, A31 * A12_ + A32 * B12_, s1, s2)
    r1 = rp(y1, 1 - Q2) * rp(u1, Q2)
    r2 = rp(y2, 1 - Q2) * rp(u2, Q2)
    g1, g2 = pat2(y1, y2, BETA1 * A21_ + BETA2 * B21_, BETA1 * A12_ + BETA2 * B12_, r1, r2)
    C21_, C12_ = prod(w1, w2)
    return pat2(y1, y2, B1 * A21_ + B2 * B21_ + B3 * C21_,
                B1 * A12_ + B2 * B12_ + B3 * C12_, g1, g2)[0]


def F(t):
    return step(Iv(t) if not isinstance(t, Iv) else t)


def dF(I):
    return step(D(I, 1)).d


def sup_abs_dF(I, pieces):
    step_ = (I.hi - I.lo) / pieces
    best = Fr(0)
    for i in range(pieces):
        d = dF(Iv(I.lo + i * step_, I.lo + (i + 1) * step_, exact=True))
        require(d.hi < 0, "F decreasing on the trap")
        best = max(best, -d.lo)
    return best


def image(I):
    return Iv(F(Iv.point(I.hi)).lo, F(Iv.point(I.lo)).hi, exact=True)


def main():
    out = certlib.out_path()
    E.reset_counters()
    params = {"scheme": "MPRK43II(gamma)", "gamma": "37/50",
              "tableau": {"a21": str(A21), "a31": str(A31), "a32": str(A32),
                          "b1": str(B1), "b2": str(B2), "b3": str(B3)},
              "patankar_exponents": {"q1": str(Q1), "q2": str(Q2)},
              "A": [["-9999000", "1000"], ["9999000", "-1000"]], "h": "1", "PREC_bits": E.PREC}
    res = {"certificate": "C7 attracting two-cycle of MPRK43II(37/50)",
           "status": "RUNNING", "parameters": params,
           "parameters_sha256": certlib.parameters_hash(params),
           "input_sha256": certlib.input_hashes("verify_mprk43_cycle.py")}
    require(min(A21, A31, A32, B1, B2, B3) > 0, "all nonzero tableau coefficients positive")
    require(A31 == Fr(73, 222) and A32 == Fr(25, 74) and Q1 == Fr(75, 74) and Q2 == Fr(3, 2), "coefficient values")
    require(B1 + B2 + B3 == 1 and A31 + A32 == A21, "consistency")
    # equilibrium
    Fp = F(Iv.point(P))
    require(Fp.contains(P), "p is a fixed point")
    dFp = dF(Iv.point(P))
    require(-1 < dFp.lo and dFp.hi < 1, "|F'(p)| < 1")
    # traps
    H = Iv(Fr(4935, 10000), Fr(5035, 10000), exact=True)
    sH = sup_abs_dF(H, 200)
    FH = image(H)
    pad = FH.width() / 20
    L = Iv(E.round_down(FH.lo - pad), E.round_up(FH.hi + pad), exact=True)
    require(L.lo > 0 and L.hi < P / 100, "L positive and far below p")
    sL = sup_abs_dF(L, 1000)
    FL = image(L)
    require(FH.inside(L) and FL.inside(H), "F swaps the traps")
    lip = sL * sH
    require(lip < 1, "contraction of F^2")
    require(H.contains(Fr(1, 2)))
    res["equilibrium"] = {"p": str(P), "F(p)": Fp.record(), "dF(p)": dFp.record()}
    res["traps"] = {"H": H.record(), "L": L.record(), "F(H)": FH.record(), "F(L)": FL.record(),
                    "sup_abs_dF_H": str(sH), "sup_abs_dF_L": str(sL),
                    "lipschitz_bound_F2": str(lip), "lipschitz_display_only": float(lip)}
    res["initial_state"] = {"y0": ["1/2", "1/2"], "in_H": True}
    res["counters"] = E.counters()
    res["floating_point_used_for_decisions"] = False
    res["status"] = "PASS"
    certlib.write_certificate(out, res)
    print("PASS  C7 MPRK43II(37/50): L=[%.4g, %.4g], H=[0.4935,0.5035], Lip(F^2) <= %.4f, F'(p) ~ %.6f"
          % (float(L.lo), float(L.hi), float(lip), float(dFp.lo)))


if __name__ == "__main__":
    main()
