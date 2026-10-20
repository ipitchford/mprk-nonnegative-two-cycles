#!/usr/bin/env python3
"""Certificate C2: complete phase portrait of the two-state witness (paper: theorems 'Counterexample' and 'Complete phase portrait').

Witness: MPRK22(4/7), A = [[-9990, 10], [9990, -10]], h = 1, i.e. K = 10^4,
p = 1/1000 in the paper's normalisation.  F is the exact scalar map of the
paper (first component on the unit-mass simplex).

The program certifies, in exact outward-rounded rational interval arithmetic:

  (M)  F'(t) < 0 for every t in (0,1).  Near t = 0 and t = 1 the derivative is
       singular (a factor t^(-1/4) or (1-t)^(-1/4)); there the program
       certifies the sign of s*F'(t) with t = s^4 (resp. 1-t = s^4), which is
       continuous up to s = 0.  On the rest of (0,1) it covers the interval
       by subintervals on which an enclosure of F' is negative.
  (S)  The sign of G(t) = F(F(t)) - t on (p, 1):
         G < 0 on (p, s_H),  G > 0 on (s_H, t_H),  G < 0 on (t_H, 1),
       with isolating rational intervals for the two zeros s_H < t_H.  Away
       from zeros the sign follows from monotonicity of F^2 and two point
       enclosures per cell; near p, s_H and t_H it follows from the sign of
       (F^2)' - 1 on an explicit interval.
  (R)  Enclosures of the cycle points t_L = F(t_H), s_L = F(s_H), of the cycle
       multipliers, of the basin length s_H - s_L, and the trapping and
       contraction facts of the counterexample theorem with sharper (subdivided) bounds.

The deductions from (M), (S), (R) to the basin theorem are written out in the
paper (proof of the phase-portrait theorem); they use only monotonicity of F^2 and the
intermediate value theorem.  No floating-point number enters a decision.
Run:  python -S -O verify_phase_portrait.py --out NEW.json
"""
from fractions import Fraction as Fr
import mprk_exact as E
from mprk_exact import Iv, require
import mprk_maps as Mm
import certlib

ALPHA = Fr(4, 7)
P = Fr(1, 1000)
KAP = Fr(1, 10000)          # kappa = 1/K
M = Mm.Method(ALPHA)
BETA, GAMMA, R = M.beta, M.gamma, M.r


def F(t):
    return Mm.F(M, t, P, KAP)


def dF(t):
    return Mm.dF(M, t, P, KAP)


def Fdec(I):
    """Image enclosure of an interval using the certified monotonicity (M)."""
    I = E.as_iv(I)
    return Iv(F(I.hi).lo, F(I.lo).hi, exact=True)


def F2(x):
    return Fdec(F(x))


# --- (M) near the endpoints -------------------------------------------------

def sdF_left(s):
    """Enclosure of s*F'(t) for t = s^4, s in an interval [0, s0] (alpha = 4/7)."""
    s = E.as_iv(s)
    t = s.ipow(4)
    a = M.alpha
    den = KAP + a
    u = (KAP * t + a * P) / den
    omu = (KAP * (1 - t) + a * (1 - P)) / den
    v1 = t / u
    v2 = (1 - t) / omu
    # v1^(3/4) = s^3 u^(-3/4),  v1^(-1/4) = u^(1/4) / s
    q1 = s.ipow(3) * u.rpow(-3, 4) * (BETA * v1 + GAMMA)
    q2 = M.Q(v2)
    dv1 = a * P / (den * u * u)
    dv2 = -(a * (1 - P)) / (den * omu * omu)
    s_dq1 = u.rpow(1, 4) * (BETA * R * v1 + GAMMA * (R - 1)) * dv1   # = s * dq1/dt
    dq2 = M.dQ(v2) * dv2
    N = KAP * t + P * q2
    D = KAP + (1 - P) * q1 + P * q2
    s_dN = s * (KAP + P * dq2)
    s_dD = (1 - P) * s_dq1 + s * (P * dq2)
    return (s_dN * D - N * s_dD) / (D * D)


def sdF_right(s):
    """Enclosure of s*F'(t) for 1 - t = s^4, s in an interval [0, s0]."""
    s = E.as_iv(s)
    w = s.ipow(4)                     # w = 1 - t
    t = 1 - w
    a = M.alpha
    den = KAP + a
    u = (KAP * t + a * P) / den
    omu = (KAP * w + a * (1 - P)) / den
    v1 = t / u
    v2 = w / omu
    q1 = M.Q(v1)
    q2 = s.ipow(3) * omu.rpow(-3, 4) * (BETA * v2 + GAMMA)
    dv1 = a * P / (den * u * u)
    dv2 = -(a * (1 - P)) / (den * omu * omu)
    dq1 = M.dQ(v1) * dv1
    s_dq2 = omu.rpow(1, 4) * (BETA * R * v2 + GAMMA * (R - 1)) * dv2
    N = KAP * t + P * q2
    D = KAP + (1 - P) * q1 + P * q2
    s_dN = s * KAP + P * s_dq2
    s_dD = s * ((1 - P) * dq1) + P * s_dq2
    return (s_dN * D - N * s_dD) / (D * D)


def F_at_one():
    """Continuous extension F(1) = kappa / (kappa + (1-p) Q(1/u(1)))."""
    a = M.alpha
    u1 = (KAP + a * P) / (KAP + a)
    q1 = M.Q(Iv.point(1) / u1)
    return KAP / (KAP + (1 - P) * q1)


def split_point(lo, hi):
    """A rational split point: geometric for wide positive ranges near 0 or 1."""
    if lo > 0 and hi / lo > 4:
        return E.round_down(Fr(E.root_bracket(lo * hi, 2)[0]))
    if hi < 1 and (1 - lo) / (1 - hi) > 4:
        g = E.root_bracket((1 - lo) * (1 - hi), 2)[0]
        return 1 - E.round_down(g)
    return (lo + hi) / 2


def cover(lo, hi, test, max_cells=200000):
    """Cover [lo, hi] by consecutive cells on which test(cell) is True."""
    cells = []
    stack = [(Fr(lo), Fr(hi))]
    while stack:
        a, b = stack.pop()
        if test(a, b):
            cells.append((a, b))
            continue
        m = split_point(a, b)
        require(a < m < b, "cover: degenerate split")
        stack.append((m, b))
        stack.append((a, m))
        require(len(cells) + len(stack) < max_cells, "cover: too many cells")
    cells.sort()
    for (a1, b1), (a2, b2) in zip(cells, cells[1:]):
        require(b1 == a2, "cover: gap")
    require(cells[0][0] == lo and cells[-1][1] == hi)
    return cells


def bisect_root(a, b, sign_at_a):
    """Shrink [a,b] (with certified opposite signs of G at the ends) to width < 2^-90 * b."""
    a, b = Fr(a), Fr(b)
    for _ in range(400):
        if b - a < b / 2 ** 90:
            break
        m = E.round_down((a + b) / 2)
        if not (a < m < b):
            m = (a + b) / 2
        g = F2(m) - m
        if g.hi < 0:
            s = -1
        elif g.lo > 0:
            s = 1
        else:
            break
        if s == sign_at_a:
            a = m
        else:
            b = m
    return a, b


def main():
    out = certlib.out_path()
    E.reset_counters()
    params = {"alpha": "4/7", "p": "1/1000", "K": "10000", "h": "1",
              "A": [["-9990", "10"], ["9990", "-10"]]}
    res = {"certificate": "C2 complete phase portrait of the two-state witness",
           "status": "RUNNING", "parameters": params,
           "parameters_sha256": certlib.parameters_hash(params),
           "input_sha256": certlib.input_hashes("verify_phase_portrait.py")}

    # ---------------------------------------------------------------- (M)
    S0 = Fr(1, 2 ** 5)                       # t <= S0^4 ~ 9.5e-7 handled by s-form
    left_cells = cover(Fr(0), S0, lambda a, b: sdF_left(Iv(a, b, exact=True)).hi < 0)
    S1 = Fr(1, 2 ** 5)                       # 1 - t <= S1^4 handled by s-form
    right_cells = cover(Fr(0), S1, lambda a, b: sdF_right(Iv(a, b, exact=True)).hi < 0)
    t_lo, t_hi = S0 ** 4, 1 - S1 ** 4
    mid_cells = cover(t_lo, t_hi, lambda a, b: dF(Iv(a, b, exact=True)).hi < 0)
    res["monotonicity"] = {
        "claim": "F'(t) < 0 for all t in (0,1)",
        "left_end": "s*F'(s^4) < 0 for s in [0, 1/32] (t in (0, 2^-20])",
        "right_end": "s*F'(1-s^4) < 0 for s in [0, 1/32] (t in [1-2^-20, 1))",
        "middle": "F'(t) < 0 on [2^-20, 1-2^-20]",
        "cells_left": len(left_cells), "cells_right": len(right_cells),
        "cells_middle": len(mid_cells)}

    # ---------------------------------------------------------------- (S)
    # neighbourhood of p: (F^2)' < 1 on [p, p + eps]
    eps = P / 20
    def deriv_F2_upper(a, b):
        I = Iv(a, b, exact=True)
        img = Fdec(I)
        d1, d2 = dF(I), dF(img)
        require(d1.hi < 0 and d2.hi < 0)
        return (d1 * d2).hi
    def deriv_F2_lower(a, b):
        I = Iv(a, b, exact=True)
        img = Fdec(I)
        d1, d2 = dF(I), dF(img)
        return (d1 * d2).lo
    near_p = cover(P, P + eps, lambda a, b: deriv_F2_upper(a, b) < 1)

    # isolating interval for s_H with (F^2)' > 1
    sH_lo, sH_hi = Fr(1080, 10000), Fr(1081, 10000)
    gl, gh = F2(sH_lo) - sH_lo, F2(sH_hi) - sH_hi
    require(gl.hi < 0 and gh.lo > 0, "sign change of G across s_H interval")
    JsH = cover(sH_lo, sH_hi, lambda a, b: deriv_F2_lower(a, b) > 1)
    # trapping interval H for t_H (counterexample theorem) and (F^2)' < 1 on it
    H_lo, H_hi = Fr(12, 25), Fr(507, 1000)
    JtH = cover(H_lo, H_hi, lambda a, b: deriv_F2_upper(a, b) < 1)
    require(F2(H_lo).lo > H_lo and F2(H_hi).hi < H_hi)

    # sign cells away from the zeros
    neg1 = cover(P + eps, sH_lo, lambda a, b: F2(b).hi < a)
    pos1 = cover(sH_hi, H_lo, lambda a, b: F2(a).lo > b)
    # last cell [X1, 1): F^2(t) <= F(F(1)) for t < 1; need F(F(1)) < X1
    F1 = F_at_one()
    require(F1.lo > 0)
    sup_F2_near_one = F(F1.lo).hi
    X1 = Fr(95, 100)
    require(sup_F2_near_one < X1, "G < 0 on [X1, 1)")
    neg2 = cover(H_hi, X1, lambda a, b: F2(b).hi < a)
    res["sign_of_G_on_(p,1)"] = {
        "near_p": {"interval": [str(P), str(P + eps)],
                   "claim": "(F^2)' < 1, hence G decreasing with G(p)=0", "cells": len(near_p)},
        "negative_cells": {"interval": [str(P + eps), str(sH_lo)], "cells": len(neg1)},
        "s_H_isolating_interval": {"interval": [str(sH_lo), str(sH_hi)],
                                   "claim": "G(lower)<0<G(upper) and (F^2)' > 1", "cells": len(JsH)},
        "positive_cells": {"interval": [str(sH_hi), str(H_lo)], "cells": len(pos1)},
        "t_H_isolating_interval": {"interval": [str(H_lo), str(H_hi)],
                                   "claim": "G(lower)>0>G(upper) and (F^2)' < 1", "cells": len(JtH)},
        "negative_cells_2": {"interval": [str(H_hi), str(X1)], "cells": len(neg2)},
        "last_cell": {"interval": [str(X1), "1"],
                      "claim": "F^2(t) < F(F(1)) < X1 <= t", "F_of_1": F1.record(),
                      "upper_bound_F2_near_1": str(sup_F2_near_one)}}

    # ---------------------------------------------------------------- (R)
    a, b = bisect_root(sH_lo, sH_hi, -1)
    sH = Iv(a, b, exact=True)
    a, b = bisect_root(H_lo, H_hi, +1)
    tH = Iv(a, b, exact=True)
    sL, tL = Fdec(sH), Fdec(tH)
    mult_att = dF(tH) * dF(tL)
    mult_rep = dF(sH) * dF(sL)
    require(0 < mult_att.lo and mult_att.hi < 1, "attracting multiplier")
    require(mult_rep.lo > 1, "repelling multiplier")
    basin = sH - sL
    # counterexample-theorem facts with monotone images and subdivided derivative bounds
    L_lo, L_hi = Fr(1, 10 ** 7), Fr(117, 10 ** 9)
    FL = Fdec(Iv(L_lo, L_hi, exact=True))
    FH = Fdec(Iv(H_lo, H_hi, exact=True))
    require(FL.inside(Iv(H_lo, H_hi, exact=True)) and FH.inside(Iv(L_lo, L_hi, exact=True)))
    n = 64
    def sup_abs_dF(lo, hi):
        step = (hi - lo) / n
        best = Fr(0)
        for i in range(n):
            d = dF(Iv(lo + i * step, lo + (i + 1) * step, exact=True))
            require(d.hi < 0)
            best = max(best, -d.lo)
        return best
    supL, supH = sup_abs_dF(L_lo, L_hi), sup_abs_dF(H_lo, H_hi)
    lip = supL * supH
    require(lip < Fr(9, 10), "contraction of F^2 on L and H below 9/10")
    res["cycles"] = {
        "attracting_cycle": {"t_H": tH.record(), "t_L": tL.record(),
                             "multiplier_F2": mult_att.record()},
        "repelling_cycle": {"s_H": sH.record(), "s_L": sL.record(),
                            "multiplier_F2": mult_rep.record()},
        "basin_of_equilibrium": {"interval": "(s_L, s_H)", "length": basin.record()},
        "local_multiplier_at_p": str(Mm.local_multiplier(ALPHA, 10000))}
    res["counterexample_traps_sharpened"] = {
        "L": [str(L_lo), str(L_hi)], "H": [str(H_lo), str(H_hi)],
        "F(L)_monotone_enclosure": FL.record(), "F(H)_monotone_enclosure": FH.record(),
        "subintervals_per_trap": n,
        "sup_abs_dF_on_L": str(supL), "sup_abs_dF_on_H": str(supH),
        "lipschitz_bound_F2": str(lip), "lipschitz_bound_F2_display_only": float(lip)}
    res["counters"] = E.counters()
    res["floating_point_used_for_decisions"] = False
    res["status"] = "PASS"
    certlib.write_certificate(out, res)
    print("PASS  C2 phase portrait;  cells (M): %d/%d/%d" % (len(left_cells), len(mid_cells), len(right_cells)))
    print("  s_H in [%.12g, %.12g], t_H in [%.12g, %.12g]" % (float(sH.lo), float(sH.hi), float(tH.lo), float(tH.hi)))
    print("  s_L ~ %.6g, t_L ~ %.6g, basin length ~ %.6g" % (float(sL.lo), float(tL.lo), float(basin.lo)))
    print("  multipliers: attracting ~ %.6g, repelling ~ %.6g; Lip(F^2) <= %.6g" %
          (float(mult_att.hi), float(mult_rep.lo), float(lip)))


if __name__ == "__main__":
    main()
