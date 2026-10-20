#!/usr/bin/env python3
"""Certificate C4: explicit parameter sets with attracting two-cycles (paper: theorem 'Explicit attracting two-cycles').

Part A.  For alpha in {13/25, 11/20, 4/7} and p in {10^-3, 10^-4, 10^-5}, the
two-state MPRK22 map F_{K,p} has an attracting positive two-cycle for EVERY
K >= 10^4, and so does the limit map K = infinity.  The step parameter enters
as kappa = 1/K in [0, 10^-4], which is covered by finitely many kappa-boxes.
Part B.  Every (alpha, p, K) in the box
    alpha in 4/7 (1 +- 10^-4),  p in 10^-3 (1 +- 10^-3),  K in 10^4 (1 +- 10^-3)
has an attracting two-cycle (an explicit form of the persistence corollary);
here alpha is an interval and powers are enclosed with log/exp series.

For each box (parameters entering as intervals) the program certifies:
    F' < 0 on L and on H,  F(H) inside int L,  F(L) inside int H,
    sup_L |F'| * sup_H |F'| < 1.
Then for every parameter in the box F^2 maps L and H into themselves and is a
contraction there, so there is a unique two-cycle with a point in each and it
attracts all of L and H; L and H are disjoint, so the period is exactly two.
The approximate cycle used to place H is computed by the program with rounded
rational arithmetic and is not part of the proof.
Run:  python -S -O verify_effective_region.py --out NEW.json
"""
from fractions import Fraction as Fr
import mprk_exact as E
from mprk_exact import Iv, require
import mprk_maps as Mm
import certlib

E.PREC = 128
ALPHA = Fr(4, 7)
M = Mm.Method(ALPHA)
P_MIN, P_MAX = Fr(1, 10 ** 5), Fr(1, 10 ** 3)
KAPPA_MAX = Fr(1, 10 ** 4)


def approx_cycle(M, p, kap):
    """Approximate attracting two-cycle (t_low, t_high) at a parameter point (not part of the proof)."""
    for start in (Fr(1, 2), Fr(9, 10), Fr(99, 100)):
        t = start
        for _ in range(600):
            t = E.round_down(Mm.F(M, t, p, kap).lo)
        t1 = E.round_down(Mm.F(M, t, p, kap).lo)
        hi, lo = max(t, t1), min(t, t1)
        if lo < p / 10 and hi > 10 * p:
            return lo, hi
    return None


def sup_abs_dF(M, I, P, K, pieces):
    lo, hi = I.lo, I.hi
    step = (hi - lo) / pieces
    best = Fr(0)
    for i in range(pieces):
        d = Mm.dF(M, Iv(lo + i * step, lo + (i + 1) * step, exact=True), P, K)
        if not d.hi < 0:
            return None
        best = max(best, -d.lo)
    return best


def image(M, I, P, K):
    """Enclosure of F(I) for all parameters in (P, K), valid when F' < 0 on I."""
    return Iv(Mm.F(M, I.hi, P, K).lo, Mm.F(M, I.lo, P, K).hi, exact=True)


def try_box(M, Mc, P, K, pc, kc):
    """Attempt a trapping-interval certificate for all parameters in the box P x K.

    H is an interval around the approximate upper cycle point; L is the
    certified image F(H) enlarged by 4% of its width on each side.  Then
    F(H) lies in int L by construction and the program checks F(L) in int H,
    F' < 0 on L and H, and sup_L|F'| sup_H|F'| < 1."""
    cyc = approx_cycle(Mc, pc, kc)
    if cyc is None:
        return None
    tl, th = cyc
    w = min(th, 1 - th)
    for rho in (Fr(1, 10), Fr(1, 30), Fr(1, 100), Fr(1, 300)):
        Hc = Iv(E.round_down(th - rho * w), E.round_up(th + rho * w), exact=True)
        dH = Mm.dF(M, Hc, P, K)
        if not dH.hi < 0:
            continue
        FH = image(M, Hc, P, K)
        pad = FH.width() / 25
        Lc = Iv(E.round_down(FH.lo - pad), E.round_up(FH.hi + pad), exact=True)
        if not (Lc.lo > 0 and Lc.hi < Hc.lo):
            continue
        sL = sup_abs_dF(M, Lc, P, K, 8)
        sH = sup_abs_dF(M, Hc, P, K, 4)
        if sL is None or sH is None or not sL * sH < 1:
            continue
        FL = image(M, Lc, P, K)
        if FL.inside(Hc) and FH.inside(Lc):
            return {"L": [str(Lc.lo), str(Lc.hi)], "H": [str(Hc.lo), str(Hc.hi)],
                    "contraction_bound": str(E.round_up(sL * sH))}, float(sL * sH)
    return None


def cover_kappa(M, p, kmax):
    """Cover kappa in [0, kmax] by certified boxes (alpha rational, p fixed)."""
    P = Iv.point(p)
    edges = [Fr(0)] + [kmax / 10 ** j for j in range(6, -1, -1)]
    stack = [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)][::-1]
    boxes, worst = [], 0.0
    while stack:
        k1, k2 = stack.pop()
        kc = k2 / 2 if k1 == 0 else E.round_down(E.root_bracket(k1 * k2, 2)[0])
        r = try_box(M, M, P, Iv(k1, k2, exact=True), p, kc)
        if r is not None:
            rec = dict(r[0]); rec["kappa"] = [str(k1), str(k2)]
            boxes.append(rec); worst = max(worst, r[1])
            continue
        require(len(boxes) + len(stack) < 100000, "too many boxes")
        km = k2 / 10 if k1 == 0 else E.round_down(E.root_bracket(k1 * k2, 2)[0])
        require(k1 < km < k2)
        stack.append((km, k2)); stack.append((k1, km))
    # the kappa boxes tile [0, kmax]
    ks = sorted((Fr(b["kappa"][0]), Fr(b["kappa"][1])) for b in boxes)
    require(ks[0][0] == 0 and ks[-1][1] == kmax and all(a[1] == b[0] for a, b in zip(ks, ks[1:])))
    return boxes, worst


def main():
    out = certlib.out_path()
    E.reset_counters()
    GRID_A = [Fr(13, 25), Fr(11, 20), Fr(4, 7)]
    GRID_P = [Fr(1, 10 ** 3), Fr(1, 10 ** 4), Fr(1, 10 ** 5)]
    params = {"part_A": {"alpha": [str(x) for x in GRID_A], "p": [str(x) for x in GRID_P],
                         "K_range": ["10000", "infinity"], "kappa_range": ["0", "1/10000"]},
              "part_B": {"alpha": ["4/7*(1-1/10000)", "4/7*(1+1/10000)"],
                         "p": ["1/1000*(1-1/1000)", "1/1000*(1+1/1000)"],
                         "K": ["10000*(1-1/1000)", "10000*(1+1/1000)"]},
              "PREC_bits": E.PREC}
    res = {"certificate": "C4 explicit parameter sets with attracting two-cycles",
           "status": "RUNNING", "parameters": params,
           "parameters_sha256": certlib.parameters_hash(params),
           "input_sha256": certlib.input_hashes("verify_effective_region.py")}
    partA = []
    nA, worstA = 0, 0.0
    for al in GRID_A:
        for pp in GRID_P:
            bx, wst = cover_kappa(Mm.Method(al), pp, KAPPA_MAX)
            partA.append({"alpha": str(al), "p": str(pp), "boxes_certified": len(bx),
                          "max_contraction_display_only": wst, "boxes": bx})
            nA += len(bx); worstA = max(worstA, wst)
            print("  alpha=%s p=%s: %d kappa-boxes, max contraction %.4f" % (al, pp, len(bx), wst), flush=True)
    res["part_A"] = {"claim": "for each listed (alpha, p): an attracting two-cycle for every K >= 10^4 and for the limit map K = infinity",
                     "cases": partA}
    Mc = Mm.Method(ALPHA)
    # part B: a three-parameter box around the witness, alpha as an interval
    a = Iv(ALPHA * (1 - Fr(1, 10 ** 4)), ALPHA * (1 + Fr(1, 10 ** 4)), exact=True)
    Mi = Mm.Method(a)
    P = Iv(Fr(1, 1000) * (1 - Fr(1, 1000)), Fr(1, 1000) * (1 + Fr(1, 1000)), exact=True)
    K = Iv(1 / (Fr(10000) * (1 + Fr(1, 1000))), 1 / (Fr(10000) * (1 - Fr(1, 1000))), exact=True)
    rB = try_box(Mi, Mc, P, K, Fr(1, 1000), Fr(1, 10000))
    require(rB is not None, "part B box")
    res["part_B"] = {"claim": "attracting two-cycle for every (alpha, p, K) in the stated box",
                     "certificate": rB[0], "contraction_display_only": rB[1]}
    res["counters"] = E.counters()
    res["floating_point_used_for_decisions"] = False
    res["status"] = "PASS"
    certlib.write_certificate(out, res)
    print("PASS  C4 part A: %d kappa-boxes in 9 cases (max contraction %.4f); part B box contraction %.4f"
          % (nA, worstA, rB[1]))


if __name__ == "__main__":
    main()
