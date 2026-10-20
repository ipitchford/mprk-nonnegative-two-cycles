"""Outward-rounded rational interval arithmetic for the MPRK22 certificates.

Python standard library only (fractions, math.isqrt is not even needed).
Every endpoint is an exact rational number.  After each arithmetic operation
the lower endpoint is rounded *down* and the upper endpoint *up* to a dyadic
rational carrying PREC significant bits, so enclosures are always valid and
numbers stay small.  Roots are computed with exact integer n-th roots and every
root enclosure is re-checked by exact rational comparison of n-th powers.
Logarithms and exponentials (used only when the method parameter alpha is
itself an interval) are enclosed by truncated series with explicit rational
remainder bounds.

All acceptance tests go through require(), which raises CertificateError and is
therefore unaffected by `python -O`.  No floating-point number is ever used
inside this module; floats appear in the verifiers only in fields labelled
"display only".
"""
from fractions import Fraction as Fr

PREC = 192          # significant bits kept after outward rounding
ROOT_CHECKS = 0     # number of exact n-th-power checks performed
SERIES_CALLS = 0    # number of log/exp series enclosures evaluated


class CertificateError(ArithmeticError):
    """Raised when an exact certificate condition fails."""


def require(condition, message="exact certificate condition failed"):
    if not condition:
        raise CertificateError(message)


def _floor_log2(x):
    """Integer e with 2**e <= |x| < 2**(e+1), for a nonzero Fraction x."""
    n, d = abs(x.numerator), x.denominator
    e = n.bit_length() - d.bit_length()
    # correct the estimate by at most one
    if e >= 0:
        if n < (d << e):
            e -= 1
    else:
        if (n << -e) < d:
            e -= 1
    return e


def round_down(x):
    x = Fr(x)
    if x == 0:
        return x
    k = PREC - _floor_log2(x)
    n, d = x.numerator, x.denominator
    if k >= 0:
        return Fr((n << k) // d, 1 << k)
    return Fr((n // (d << -k)) << -k)


def round_up(x):
    return -round_down(-Fr(x))


def iroot(N, n):
    """Floor of the real n-th root of the nonnegative integer N."""
    require(N >= 0 and n >= 1)
    if N < 2 or n == 1:
        return N
    x = 1 << ((N.bit_length() + n - 1) // n)
    while True:
        y = ((n - 1) * x + N // x ** (n - 1)) // n
        if y >= x:
            break
        x = y
    while x ** n > N:
        x -= 1
    while (x + 1) ** n <= N:
        x += 1
    return x


def root_bracket(x, n):
    """Rationals lo <= x**(1/n) < hi, checked exactly; x >= 0 rational."""
    global ROOT_CHECKS
    x = Fr(x)
    require(x >= 0 and n >= 1)
    if x == 0:
        return Fr(0), Fr(0)
    m = PREC - _floor_log2(x) // n + 2
    if m >= 0:
        N = (x.numerator << (n * m)) // x.denominator
        k = iroot(N, n)
        lo, hi = Fr(k, 1 << m), Fr(k + 1, 1 << m)
    else:
        N = x.numerator // (x.denominator << (-n * m))
        k = iroot(N, n)
        lo, hi = Fr(k << -m), Fr((k + 1) << -m)
    require(lo ** n <= x < hi ** n, "n-th root bracket")
    ROOT_CHECKS += 1
    return lo, hi


class Iv:
    """Closed interval [lo, hi] with exact rational endpoints."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None, exact=False):
        lo = Fr(lo)
        hi = lo if hi is None else Fr(hi)
        require(lo <= hi, "interval endpoints out of order")
        if exact:
            self.lo, self.hi = lo, hi
        else:
            self.lo, self.hi = round_down(lo), round_up(hi)

    # --- construction helpers -------------------------------------------
    @staticmethod
    def point(x):
        return Iv(x, x, exact=True)

    # --- predicates ---------------------------------------------------------
    def contains(self, x):
        return self.lo <= x <= self.hi

    def inside(self, other):
        """Strictly inside the interior of `other`."""
        other = as_iv(other)
        return other.lo < self.lo and self.hi < other.hi

    def subset(self, other):
        other = as_iv(other)
        return other.lo <= self.lo and self.hi <= other.hi

    def positive(self):
        return self.lo > 0

    def negative(self):
        return self.hi < 0

    # --- arithmetic -----------------------------------------------------------
    def __add__(self, o):
        o = as_iv(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return Iv(-self.hi, -self.lo, exact=True)

    def __sub__(self, o):
        return self + (-as_iv(o))

    def __rsub__(self, o):
        return as_iv(o) + (-self)

    def __mul__(self, o):
        o = as_iv(o)
        c = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Iv(min(c), max(c))

    __rmul__ = __mul__

    def recip(self):
        require(self.lo > 0 or self.hi < 0, "division by an interval containing 0")
        return Iv(1 / self.hi, 1 / self.lo)

    def __truediv__(self, o):
        return self * as_iv(o).recip()

    def __rtruediv__(self, o):
        return as_iv(o) * self.recip()

    def ipow(self, k):
        """Integer power k >= 0 (exact monotone rule)."""
        require(isinstance(k, int) and k >= 0)
        if k == 0:
            return Iv.point(1)
        if self.lo >= 0:
            return Iv(self.lo ** k, self.hi ** k)
        if self.hi <= 0:
            a, b = (-self.hi) ** k, (-self.lo) ** k
            return Iv(a, b) if k % 2 == 0 else Iv(-b, -a)
        m = max(-self.lo, self.hi) ** k
        return Iv(0, m) if k % 2 == 0 else Iv(self.lo ** k, self.hi ** k)

    def root(self, n):
        """Principal n-th root of a nonnegative interval."""
        require(self.lo >= 0, "root of an interval with negative part")
        return Iv(root_bracket(self.lo, n)[0], root_bracket(self.hi, n)[1])

    def rpow(self, k, m):
        """self ** (k/m) for a positive interval and integers k, m > 0 (k may be negative)."""
        require(self.lo > 0 and m > 0)
        if k >= 0:
            return self.ipow(k).root(m)
        return self.ipow(-k).root(m).recip()

    def hull(self, o):
        o = as_iv(o)
        return Iv(min(self.lo, o.lo), max(self.hi, o.hi), exact=True)

    def mid(self):
        return (self.lo + self.hi) / 2

    def width(self):
        return self.hi - self.lo

    def record(self):
        return {"lower": str(self.lo), "upper": str(self.hi),
                "decimal_lower_display_only": float(self.lo),
                "decimal_upper_display_only": float(self.hi)}

    def __repr__(self):
        return "Iv(%.6g, %.6g)" % (float(self.lo), float(self.hi))


def as_iv(x):
    return x if isinstance(x, Iv) else Iv.point(x)


# ---------------------------------------------------------------------------
# logarithm and exponential (for interval-valued alpha)
# ---------------------------------------------------------------------------
_LN2 = None


def _atanh_series(s, terms):
    """Enclosure of 2*atanh(s) for a rational 0 <= s <= 1/5 with remainder."""
    S = Iv.point(s)
    S2 = S * S
    total = Iv.point(0)
    power = S
    for k in range(terms):
        total = total + power / (2 * k + 1)
        power = power * S2
    # remainder of sum_{k>=terms} s^(2k+1)/(2k+1) <= s^(2terms+1)/((2terms+1)(1-s^2))
    rem = Fr(s) ** (2 * terms + 1) / ((2 * terms + 1) * (1 - Fr(s) ** 2))
    return 2 * (total + Iv(0, rem))


def ln2():
    global _LN2
    if _LN2 is None:
        # ln 2 = 2 atanh(1/3); use s = 1/3 with enough terms (each term ~3.17 bits)
        S = Iv.point(Fr(1, 3))
        S2 = S * S
        total = Iv.point(0)
        power = S
        terms = PREC // 3 + 8
        for k in range(terms):
            total = total + power / (2 * k + 1)
            power = power * S2
        s = Fr(1, 3)
        rem = s ** (2 * terms + 1) / ((2 * terms + 1) * (1 - s * s))
        _LN2 = 2 * (total + Iv(0, rem))
    return _LN2


def ln_point(x):
    """Enclosure of log(x) for a rational x > 0."""
    global SERIES_CALLS
    x = Fr(x)
    require(x > 0)
    SERIES_CALLS += 1
    e = _floor_log2(x)
    m = x / Fr(2) ** e            # 1 <= m < 2
    if m * m > 2:                 # move m into [1/sqrt2, sqrt2]
        m /= 2
        e += 1
    s = (m - 1) / (m + 1)         # |s| <= 0.1716
    s = abs(s)
    terms = PREC // 5 + 6
    val = _atanh_series(s, terms)
    if m < 1:
        val = -val
    return val + e * ln2()


def exp_point(y):
    """Enclosure of exp(y) for a rational y."""
    global SERIES_CALLS
    y = Fr(y)
    SERIES_CALLS += 1
    k = (y * Fr(14427, 10000) + Fr(1, 2)).__floor__()   # any integer is valid
    f = Iv.point(y) - k * ln2()                          # |f| < 0.36
    bound = max(abs(f.lo), abs(f.hi))
    require(bound < Fr(1, 2))
    total = Iv.point(1)
    term = Iv.point(1)
    n = 1
    while True:
        term = term * f / n
        total = total + term
        n += 1
        # remainder after n-1 terms: |f|^n/n! * e^{|f|} <= 2 |f|^n / n!
        tail = Fr(2) * bound ** n
        fact = 1
        for j in range(2, n + 1):
            fact *= j
        tail /= fact
        if tail < Fr(1, 2 ** (PREC + 8)):
            break
    total = total + Iv(-tail, tail)
    scale = Fr(2) ** k
    return Iv(total.lo * scale, total.hi * scale)


def ln_iv(v):
    v = as_iv(v)
    require(v.lo > 0)
    if v.lo == v.hi:
        return ln_point(v.lo)
    return Iv(ln_point(v.lo).lo, ln_point(v.hi).hi, exact=True)


def exp_iv(v):
    v = as_iv(v)
    if v.lo == v.hi:
        return exp_point(v.lo)
    return Iv(exp_point(v.lo).lo, exp_point(v.hi).hi, exact=True)


def pow_iv(v, e):
    """v ** e for a positive interval v and an interval exponent e."""
    return exp_iv(as_iv(e) * ln_iv(v))


def counters():
    return {"exact_root_checks": ROOT_CHECKS, "log_exp_series_enclosures": SERIES_CALLS}


def reset_counters():
    global ROOT_CHECKS, SERIES_CALLS
    ROOT_CHECKS = 0
    SERIES_CALLS = 0
