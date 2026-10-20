"""Exact interval enclosures of the MPRK22(alpha) maps used by the certificates.

Standard library only; built on mprk_exact.  Notation follows the paper.

Two-state scalar map (Section 2 of the paper).  The test matrix is
K*A_p with A_p = [[-(1-p), p], [1-p, -p]], and the step-size parameter is
written kappa = 1/K so that the limit K -> infinity is the regular value
kappa = 0:

    u(t)   = (kappa*t + alpha*p) / (kappa + alpha)        (stage, first component)
    1-u(t) = (kappa*(1-t) + alpha*(1-p)) / (kappa + alpha)
    v1 = t/u,   v2 = (1-t)/(1-u)
    Q(v) = v^(r-1) * (beta*v + gamma),  r = 1/alpha, beta = 1-1/(2alpha), gamma = 1/(2alpha)
    F(t) = (kappa*t + p*Q(v2)) / (kappa + (1-p)*Q(v1) + p*Q(v2)).

Multiplying numerator and denominator by K gives the paper's formula
F = (t + K p q2)/(1 + K((1-p) q1 + p q2)) with q1 = Q(v1), q2 = Q(v2).

alpha may be a Fraction m/n (exact rational powers via integer roots) or an
Iv (interval exponent via log/exp enclosures).
"""
from fractions import Fraction as Fr
from mprk_exact import Iv, as_iv, require, pow_iv


class Method:
    """Coefficients of MPRK22(alpha) for rational or interval alpha."""

    def __init__(self, alpha):
        if isinstance(alpha, Iv):
            require(alpha.lo >= Fr(1, 2))
            self.rational = False
            self.alpha = alpha
            self.r = alpha.recip()
            self.gamma = Iv.point(1) / (2 * alpha)
            self.beta = 1 - self.gamma
        else:
            alpha = Fr(alpha)
            require(alpha >= Fr(1, 2))
            self.rational = True
            self.alpha = Iv.point(alpha)
            self.m, self.n = alpha.numerator, alpha.denominator   # alpha = m/n, r = n/m
            self.r = Iv.point(Fr(self.n, self.m))
            self.gamma = Iv.point(Fr(self.n, 2 * self.m))
            self.beta = Iv.point(1 - Fr(self.n, 2 * self.m))

    # v^(r-1) and v^(r-2) for positive intervals v
    def pow_rm1(self, v):
        if self.rational:
            return v.rpow(self.n - self.m, self.m)
        return pow_iv(v, self.r - 1)

    def pow_rm2(self, v):
        if self.rational:
            return v.rpow(self.n - 2 * self.m, self.m)
        return pow_iv(v, self.r - 2)

    def Q(self, v):
        return self.pow_rm1(v) * (self.beta * v + self.gamma)

    def dQ(self, v):
        """Q'(v) = v^(r-2) (beta r v + gamma (r-1))."""
        return self.pow_rm2(v) * (self.beta * self.r * v + self.gamma * (self.r - 1))


def _stage(M, t, p, kap):
    a = M.alpha
    den = kap + a
    u = (kap * t + a * p) / den
    omu = (kap * (1 - t) + a * (1 - p)) / den
    return u, omu, den


def F(M, t, p, kap):
    """Enclosure of the scalar MPRK22 map for t, p, kap intervals (or rationals)."""
    t, p, kap = as_iv(t), as_iv(p), as_iv(kap)
    require(t.lo > 0 and t.hi < 1 and p.lo > 0 and p.hi < 1 and kap.lo >= 0)
    u, omu, _ = _stage(M, t, p, kap)
    q1 = M.Q(t / u)
    q2 = M.Q((1 - t) / omu)
    return (kap * t + p * q2) / (kap + (1 - p) * q1 + p * q2)


def dF(M, t, p, kap):
    """Enclosure of dF/dt."""
    t, p, kap = as_iv(t), as_iv(p), as_iv(kap)
    require(t.lo > 0 and t.hi < 1 and p.lo > 0 and p.hi < 1 and kap.lo >= 0)
    a = M.alpha
    u, omu, den = _stage(M, t, p, kap)
    v1 = t / u
    v2 = (1 - t) / omu
    q1, q2 = M.Q(v1), M.Q(v2)
    dv1 = a * p / (den * u * u)
    dv2 = -(a * (1 - p)) / (den * omu * omu)
    dq1 = M.dQ(v1) * dv1
    dq2 = M.dQ(v2) * dv2
    N = kap * t + p * q2
    D = kap + (1 - p) * q1 + p * q2
    dN = kap + p * dq2
    dD = (1 - p) * dq1 + p * dq2
    return (dN * D - N * dD) / (D * D)


def local_multiplier(alpha, K):
    """Exact F'(p) = R(-K) = (-K^2 + 2 alpha K + 2) / (2 (1 + alpha K)(1 + K))."""
    alpha, K = Fr(alpha), Fr(K)
    return (-K * K + 2 * alpha * K + 2) / (2 * (1 + alpha * K) * (1 + K))


# ---------------------------------------------------------------------------
# N-state linear map for rational alpha (used for the three-state example)
# ---------------------------------------------------------------------------

def solve(Mtx, b):
    """Interval Gaussian elimination without pivoting (used on M-matrices).

    Mtx is a list of lists of Iv, b a list of Iv.  The diagonal pivots are
    required to stay positive, which holds for the nonsingular M-matrices that
    occur here; the enclosure is valid whenever the pivots exclude zero.
    """
    n = len(b)
    A = [row[:] for row in Mtx]
    x = b[:]
    for k in range(n):
        require(A[k][k].lo > 0, "pivot must be positive")
        for i in range(k + 1, n):
            f = A[i][k] / A[k][k]
            for j in range(k + 1, n):
                A[i][j] = A[i][j] - f * A[k][j]
            x[i] = x[i] - f * x[k]
    sol = [None] * n
    for i in reversed(range(n)):
        s = x[i]
        for j in range(i + 1, n):
            s = s - A[i][j] * sol[j]
        sol[i] = s / A[i][i]
    return sol


class LinearSystem:
    """MPRK22(alpha) on y' = A y for a rational conservative Metzler matrix A."""

    def __init__(self, alpha, A, h):
        self.M = Method(Fr(alpha))
        self.N = len(A)
        self.A = [[Fr(x) for x in row] for row in A]
        self.h = Fr(h)
        for j in range(self.N):
            require(sum(self.A[i][j] for i in range(self.N)) == 0, "column sums must vanish")
            for i in range(self.N):
                if i != j:
                    require(self.A[i][j] >= 0, "off-diagonal entries must be nonnegative")
        a = self.M.alpha.lo
        # stage matrix I - alpha h A is exact; its inverse P is computed exactly
        S = [[(1 if i == j else 0) - a * self.h * self.A[i][j] for j in range(self.N)]
             for i in range(self.N)]
        self.P = exact_inverse(S)

    def stage(self, y):
        return [sum(self.P[i][j] * y[j] for j in range(self.N)) for i in range(self.N)]

    def weights(self, y, u):
        return [self.M.Q(y[i] / u[i]) for i in range(self.N)]

    def step(self, y):
        """Interval enclosure of one MPRK22 step for an interval vector y > 0."""
        y = [as_iv(v) for v in y]
        u = self.stage(y)
        q = self.weights(y, u)
        hA = self.h
        Mtx = [[(Iv.point(1) if i == j else Iv.point(0)) - hA * self.A[i][j] * q[j]
                for j in range(self.N)] for i in range(self.N)]
        return solve(Mtx, y)

    def jacobian(self, y):
        """Interval enclosure of the full N x N derivative of one step at y.

        With M(q) = I - h A diag(q), y+ = M^{-1} y, u = P y and
        q_i = Q(y_i/u_i):  dy+ = M^{-1} (dy + h A diag(y+) dq),
        dq_i = Q'(v_i) (dy_i/u_i - y_i du_i/u_i^2),  du = P dy.
        """
        N = self.N
        y = [as_iv(v) for v in y]
        u = self.stage(y)
        v = [y[i] / u[i] for i in range(N)]
        q = [self.M.Q(v[i]) for i in range(N)]
        dQ = [self.M.dQ(v[i]) for i in range(N)]
        Mtx = [[(Iv.point(1) if i == j else Iv.point(0)) - self.h * self.A[i][j] * q[j]
                for j in range(N)] for i in range(N)]
        yp = solve(Mtx, y)
        # dq/dy matrix
        Dq = [[None] * N for _ in range(N)]
        for i in range(N):
            for j in range(N):
                term = -(y[i] * self.P[i][j]) / (u[i] * u[i])
                if i == j:
                    term = term + Iv.point(1) / u[i]
                Dq[i][j] = dQ[i] * term
        # right-hand sides: column j of  I + h A diag(y+) Dq
        J = [[None] * N for _ in range(N)]
        for j in range(N):
            rhs = []
            for i in range(N):
                s = Iv.point(1 if i == j else 0)
                for k in range(N):
                    if self.A[i][k] != 0:
                        s = s + self.h * self.A[i][k] * yp[k] * Dq[k][j]
                rhs.append(s)
            col = solve(Mtx, rhs)
            for i in range(N):
                J[i][j] = col[i]
        return J


def exact_inverse(S):
    """Exact rational inverse by Gauss-Jordan elimination."""
    n = len(S)
    A = [[Fr(x) for x in row] + [Fr(1 if i == j else 0) for j in range(n)] for i, row in enumerate(S)]
    for k in range(n):
        piv = next(i for i in range(k, n) if A[i][k] != 0)
        A[k], A[piv] = A[piv], A[k]
        pk = A[k][k]
        A[k] = [x / pk for x in A[k]]
        for i in range(n):
            if i != k and A[i][k] != 0:
                f = A[i][k]
                A[i] = [x - f * z for x, z in zip(A[i], A[k])]
    return [row[n:] for row in A]
