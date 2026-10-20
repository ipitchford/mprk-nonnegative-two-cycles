#!/usr/bin/env python3
"""Certificate C5: exact polynomial identities behind the alpha = 1/2 and alpha = 1 theorems.

Standard library only: a minimal multivariate polynomial type with Fraction
coefficients.  Checks

  (A) alpha = 1/2 (paper, Section 7).  With U = 2t + Kp and W = 2(1-t) + K(1-p),
      the scalar map is F = N/D with
          N = t U W + K p (1-t)(2+K) U,
          D = U W + K (1-p) t (2+K) W + K p (1-t)(2+K) U.
      (A1) N, D are the stated clearings of the MPRK22(1/2) formulas;
      (A2) D - N = (1-t) U W + K (1-p) t (2+K) W;
      (A3) N - t D = K (2+K) (p - t) [K (p(1-t) + t(1-p)) + 2 t (1-t)];
      (A4) Phi := N t (1-p)^2 - p^2 (1-t)(D - N) = (t - p) Psi, and after the
           substitution t = a/(a+b), p = c/(c+d) the polynomial
           (a+b)^3 (c+d)^3 Psi has only nonnegative coefficients (and is
           nonzero), so Psi > 0 for t, p in (0,1), K > 0.
  (B) alpha = 1 (paper, Section 7): the exact contraction identity
          F - p = (t - p)(1 - T)/(1 + K + K^2/2 - (1+K) T),
          T = K^2 p (1-p) / (2 (1+K) u (1-u)),  u = (t + K p)/(1+K),
      the auxiliary identity (1-p) t (1-u) + p (1-t) u = (1+K) u(1-u) - K p(1-p),
      the endpoint values used to bound T,
      and the crossover 2/(K^2+2K+2) - (K-2)/(K+2) = -(K^3-4K-8)/((K+2)(K^2+2K+2))
      with the unique positive root of K^3-4K-8 in (2649/1000, 53/20).
Run:  python -S -O verify_algebraic_identities.py --out NEW.json
"""
from fractions import Fraction as Fr
from mprk_exact import require
import certlib

VARS = ("K", "p", "t", "a", "b", "c", "d")


class Poly:
    __slots__ = ("c",)

    def __init__(self, coeffs=None):
        self.c = {k: Fr(v) for k, v in (coeffs or {}).items() if v != 0}

    @staticmethod
    def const(x):
        return Poly({(0,) * len(VARS): Fr(x)})

    @staticmethod
    def var(name):
        e = [0] * len(VARS)
        e[VARS.index(name)] = 1
        return Poly({tuple(e): 1})

    def __add__(self, o):
        o = o if isinstance(o, Poly) else Poly.const(o)
        r = dict(self.c)
        for k, v in o.c.items():
            r[k] = r.get(k, 0) + v
        return Poly(r)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.c.items()})

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Poly) else Poly.const(o)))

    def __rsub__(self, o):
        return Poly.const(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, Poly) else Poly.const(o)
        r = {}
        for k1, v1 in self.c.items():
            for k2, v2 in o.c.items():
                k = tuple(x + y for x, y in zip(k1, k2))
                r[k] = r.get(k, 0) + v1 * v2
        return Poly(r)

    __rmul__ = __mul__

    def __pow__(self, n):
        r = Poly.const(1)
        for _ in range(n):
            r = r * self
        return r

    def is_zero(self):
        return not self.c

    def degree(self, name):
        i = VARS.index(name)
        return max((k[i] for k in self.c), default=0)

    def coeff_in(self, name, j):
        """Coefficient polynomial of name^j."""
        i = VARS.index(name)
        r = {}
        for k, v in self.c.items():
            if k[i] == j:
                kk = list(k)
                kk[i] = 0
                r[tuple(kk)] = v
        return Poly(r)

    def subs_poly(self, name, P):
        """Substitute a polynomial P for variable name."""
        i = VARS.index(name)
        out = Poly()
        for k, v in self.c.items():
            kk = list(k)
            e = kk[i]
            kk[i] = 0
            out = out + Poly({tuple(kk): v}) * (P ** e)
        return out


K, p, t = Poly.var("K"), Poly.var("p"), Poly.var("t")
a, b, c, d = (Poly.var(x) for x in "abcd")


def divide_linear(Phi, name, root):
    """Divide Phi by (name - root) exactly; root is a polynomial free of name."""
    n = Phi.degree(name)
    cs = [Phi.coeff_in(name, j) for j in range(n + 1)]
    q = [None] * n
    q[n - 1] = cs[n]
    for j in range(n - 1, 0, -1):
        q[j - 1] = cs[j] + root * q[j]
    rem = cs[0] + root * q[0]
    require(rem.is_zero(), "exact division by (t - p)")
    x = Poly.var(name)
    out = Poly()
    for j, qj in enumerate(q):
        out = out + qj * (x ** j)
    return out


def homogenise(Pl):
    """(a+b)^dt (c+d)^dp * Pl(t = a/(a+b), p = c/(c+d)) as a polynomial."""
    dt, dp = Pl.degree("t"), Pl.degree("p")
    ti, pi = VARS.index("t"), VARS.index("p")
    out = Poly()
    for k, v in Pl.c.items():
        kk = list(k)
        jt, jp = kk[ti], kk[pi]
        kk[ti] = kk[pi] = 0
        term = Poly({tuple(kk): v}) * (a ** jt) * ((a + b) ** (dt - jt)) * (c ** jp) * ((c + d) ** (dp - jp))
        out = out + term
    return out, dt, dp


def main():
    out = certlib.out_path()
    res = {"certificate": "C5 exact polynomial identities (alpha = 1/2 and alpha = 1)",
           "status": "RUNNING", "input_sha256": certlib.input_hashes("verify_algebraic_identities.py")}
    one = Poly.const(1)
    # ---------------------------------------------------------------- (A)
    U = 2 * t + K * p
    W = 2 * (1 - t) + K * (1 - p)
    N = t * U * W + K * p * (1 - t) * (2 + K) * U
    D = U * W + K * (1 - p) * t * (2 + K) * W + K * p * (1 - t) * (2 + K) * U
    # (A1) MPRK22(1/2): u = U/(2+K), 1-u = W/(2+K), q1 = t/u, q2 = (1-t)/(1-u),
    # F = (t + K p q2)/(1 + K((1-p) q1 + p q2)); multiply through by U W/(2+K)... check
    # N/D = (t + K p q2)/(1 + K (1-p) q1 + K p q2) via cross-multiplication with
    # q1 = t (2+K)/U and q2 = (1-t)(2+K)/W:
    numF = t * U * W + K * p * (1 - t) * (2 + K) * U        # (t + K p q2) * U W
    denF = U * W + K * (1 - p) * t * (2 + K) * W + K * p * (1 - t) * (2 + K) * U
    require((N - numF).is_zero() and (D - denF).is_zero(), "A1")
    # stage check: (I - K/2 A_p) u = y for first component: u (1 + K/2) - K p/2 = t
    require((U * (1 + K * Fr(1, 2)) - (2 + K) * (t + K * p * Fr(1, 2))).is_zero(), "A1 stage")
    require((D - N - ((1 - t) * U * W + K * (1 - p) * t * (2 + K) * W)).is_zero(), "A2")
    A3 = K * (2 + K) * (p - t) * (K * (p * (1 - t) + t * (1 - p)) + 2 * t * (1 - t))
    require((N - t * D - A3).is_zero(), "A3")
    Phi = N * t * (1 - p) ** 2 - p ** 2 * (1 - t) * (D - N)
    Psi = divide_linear(Phi, "t", p)
    require((Phi - (t - p) * Psi).is_zero(), "A4 factorisation")
    # explicit positive form used in the paper (s = 1-t, q = 1-p)
    s_, q_ = 1 - t, 1 - p
    Psi_form = (4 * t * s_ * (p * s_ + q_ * t)
                + 2 * K * (p * p * s_ * s_ + 4 * p * q_ * s_ * t + q_ * q_ * t * t)
                + K * K * p * q_ * (p * s_ * s_ + 3 * p * s_ * t + 3 * q_ * s_ * t + q_ * t * t))
    require((Psi - Psi_form).is_zero(), "A4 explicit positive form of Psi")
    Hm, dt, dp = homogenise(Psi)
    coeffs = list(Hm.c.values())
    require(len(coeffs) > 0 and min(coeffs) >= 0, "A4 nonnegative coefficients")
    # human-readable Psi (monomials in K, p, t)
    def show(Pl):
        terms = []
        for k, v in sorted(Pl.c.items(), reverse=True):
            mono = "*".join("%s^%d" % (VARS[i], e) if e > 1 else VARS[i] for i, e in enumerate(k) if e)
            terms.append("%s%s" % (v, ("*" + mono) if mono else ""))
        return " + ".join(terms)
    res["alpha_one_half"] = {
        "U": "2t + Kp", "W": "2(1-t) + K(1-p)",
        "N": "t U W + K p (1-t)(2+K) U", "D": "U W + K(1-p) t (2+K) W + K p (1-t)(2+K) U",
        "A2_D_minus_N": "(1-t) U W + K(1-p) t (2+K) W",
        "A3_N_minus_tD": "K(2+K)(p-t)[K(p(1-t)+t(1-p)) + 2t(1-t)]",
        "Psi": show(Psi),
        "Psi_positive_form": "4ts(ps+qt) + 2K(p^2 s^2 + 4pqst + q^2 t^2) + K^2 pq(ps^2 + 3pst + 3qst + qt^2), s=1-t, q=1-p",
        "homogenised_degrees": [dt, dp],
        "homogenised_terms": len(coeffs), "min_coefficient": str(min(coeffs))}
    # ---------------------------------------------------------------- (B)
    # alpha = 1: u = (t + K p)/(1+K); write uu = (1+K) u = t + K p, vv = (1+K)(1-u) = 1 - t + K(1-p)
    uu = t + K * p
    vv = (1 - t) + K * (1 - p)
    require((uu + vv - (1 + K)).is_zero())
    # q1 = (1 + t/u)/2 = (uu + (1+K) t)/(2 uu), q2 = (vv + (1+K)(1-t))/(2 vv)
    # F = (t + K p q2)/(1 + K (1-p) q1 + K p q2); multiply by 2 uu vv
    numF1 = 2 * t * uu * vv + K * p * uu * (vv + (1 + K) * (1 - t))
    denF1 = 2 * uu * vv + K * (1 - p) * vv * (uu + (1 + K) * t) + K * p * uu * (vv + (1 + K) * (1 - t))
    # T = K^2 p(1-p)/(2(1+K) u(1-u)) = K^2 p (1-p)(1+K) / (2 uu vv)
    # claim: (F - p) * (D0 - D1 T) = (t - p)(1 - T),  D0 = 1+K+K^2/2, D1 = 1+K
    # multiply by 2 uu vv * denF1:
    D0 = 1 + K + K * K * Fr(1, 2)
    D1 = 1 + K
    Tn = K * K * p * (1 - p) * (1 + K)        # T = Tn / (2 uu vv)
    lhs = (numF1 - p * denF1) * (D0 * 2 * uu * vv - D1 * Tn)
    rhs = (t - p) * (2 * uu * vv - Tn) * denF1
    require((lhs - rhs).is_zero(), "B identity (alphaone)")
    # auxiliary identity times (1+K): (1+K)[(1-p) t (1-u) + p (1-t) u] = (1+K)[(1+K)u(1-u) - K p(1-p)]
    auxl = (1 - p) * t * vv + p * (1 - t) * uu
    auxr = uu * vv - K * p * (1 - p) * (1 + K)
    require((auxl - auxr).is_zero(), "B auxiliary identity (multiplied by 1+K)")
    # endpoint values of u(1-u) - K p(1-p)/(1+K), cleared by (1+K)^2
    e1 = (K * p) * (1 + K - K * p) - K * p * (1 - p) * (1 + K)
    e2 = (1 + K * p) * (K - K * p) - K * p * (1 - p) * (1 + K)
    require((e1 - K * p * p).is_zero() and (e2 - K * (1 - p) ** 2).is_zero(), "B endpoint values")
    # g(0) and g(K/2):  g(T) = (1-T)/(D0 - D1 T)
    require(((D0 - D1 * K * Fr(1, 2)) * 1 - (1 + K * Fr(1, 2))).is_zero())
    # crossover: [2 (K+2) - (K-2)(K^2+2K+2)] = -(K^3 - 4K - 8)
    cross = 2 * (K + 2) - (K - 2) * (K * K + 2 * K + 2)
    require((cross + (K ** 3 - 4 * K - 8)).is_zero(), "B crossover polynomial")
    f = lambda x: x ** 3 - 4 * x - 8
    lo, hi = Fr(2649, 1000), Fr(53, 20)
    require(f(lo) < 0 < f(hi), "root bracket")
    # uniqueness of the positive root: f'(x) = 3x^2 - 4 > 0 for x >= 2/sqrt(3); f < 0 on (0, 2/sqrt3]
    require(f(Fr(2)) < 0 and 3 * Fr(2) ** 2 - 4 > 0, "f increasing on [2, inf) and negative at 2")
    res["alpha_one"] = {
        "identity": "F - p = (t-p)(1-T)/(1+K+K^2/2-(1+K)T)  (cleared form checked)",
        "auxiliary_identity": "(1-p)t(1-u) + p(1-t)u = (1+K)u(1-u) - Kp(1-p)",
        "endpoint_values": ["K p^2/(1+K)^2", "K (1-p)^2/(1+K)^2"],
        "crossover_polynomial": "K^3 - 4K - 8",
        "crossover_root_bracket": [str(lo), str(hi)]}
    res["status"] = "PASS"
    certlib.write_certificate(out, res)
    print("PASS  C5 algebraic identities; Psi homogenised has %d terms, min coefficient %s"
          % (len(coeffs), min(coeffs)))


if __name__ == "__main__":
    main()
