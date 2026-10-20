#!/usr/bin/env python3
"""Certificate C6: a three-parameter region in which global convergence fails (paper: theorem 'An explicit region of non-convergence').

Claim.  For every alpha in [53/100, 56/100], every p in [10^-4, 10^-3] and
every K in [10^4, 10^6], the two-state MPRK22 map F_{K,p} has a positive
periodic orbit of prime period two.  Consequently MPRK22(alpha) is not
globally convergent for these parameters.

Argument (written out in the paper).  For finite K the map F extends
continuously to [0,1] with F(1) in (0,1), so G(t) = F(F(t)) - t satisfies
G(1) < 0.  If G(x) > 0 for some rational x > p, the intermediate value theorem
gives a zero of G in (x,1); it is not the unique fixed point p (Lemma 1),
hence it lies on a two-cycle.  The program covers the parameter box by
sub-boxes (alpha, p, K all as intervals; real powers via log/exp series with
rational remainder bounds) and certifies, on each sub-box, G(x) > 0 for one
test point x = 2^-j with x > 2 p_max.
Run:  python -S -O verify_existence_region.py --out NEW.json
"""
from fractions import Fraction as Fr
import mprk_exact as E
from mprk_exact import Iv, require
import mprk_maps as Mm
import certlib

E.PREC = 96
A_MIN, A_MAX = Fr(53, 100), Fr(56, 100)
P_MIN, P_MAX = Fr(1, 10 ** 4), Fr(1, 10 ** 3)
K_MIN, K_MAX = Fr(10 ** 4), Fr(10 ** 6)


def test_box(a1, a2, p1, p2, K1, K2):
    M = Mm.Method(Iv(a1, a2, exact=True))
    P = Iv(p1, p2, exact=True)
    kap = Iv(1 / K2, 1 / K1, exact=True)
    j = 1
    while Fr(1, 2 ** j) > 2 * p2:
        x = Fr(1, 2 ** j)
        y = Mm.F(M, x, P, kap)
        if y.lo > 0 and y.hi < 1:
            z = Mm.F(M, y, P, kap)
            if z.lo > x:
                return j
        j += 1
    return None


def main():
    out = certlib.out_path()
    E.reset_counters()
    params = {"alpha": [str(A_MIN), str(A_MAX)], "p": [str(P_MIN), str(P_MAX)],
              "K": [str(K_MIN), str(K_MAX)], "PREC_bits": E.PREC}
    res = {"certificate": "C6 region of parameters with a two-cycle (non-convergence)",
           "status": "RUNNING", "parameters": params,
           "parameters_sha256": certlib.parameters_hash(params),
           "input_sha256": certlib.input_hashes("verify_existence_region.py")}
    a_edges = [A_MIN + i * Fr(1, 400) for i in range(int((A_MAX - A_MIN) / Fr(1, 400)) + 1)]
    require(a_edges[-1] == A_MAX)
    p_edges = [P_MIN]
    while p_edges[-1] < P_MAX:
        p_edges.append(min(P_MAX, p_edges[-1] * Fr(5, 4)))
    k_edges = [K_MIN]
    while k_edges[-1] < K_MAX:
        k_edges.append(min(K_MAX, k_edges[-1] * Fr(5, 4)))
    stack = [(a_edges[i], a_edges[i + 1], p_edges[j], p_edges[j + 1], k_edges[k], k_edges[k + 1])
             for i in range(len(a_edges) - 1) for j in range(len(p_edges) - 1)
             for k in range(len(k_edges) - 1)][::-1]
    boxes = []
    while stack:
        a1, a2, p1, p2, K1, K2 = stack.pop()
        j = test_box(a1, a2, p1, p2, K1, K2)
        if j is not None:
            boxes.append([str(a1), str(a2), str(p1), str(p2), str(K1), str(K2), j])
            if len(boxes) % 500 == 0:
                print("  %d sub-boxes certified, %d pending" % (len(boxes), len(stack)), flush=True)
            continue
        require(len(boxes) + len(stack) < 100000, "too many boxes")
        # split the relatively widest coordinate (alpha width in units of 1/400)
        wa = (a2 - a1) * 1600
        wp = Fr(p2, p1)
        wk = Fr(K2, K1)
        if wa >= 1 and wa * wa >= max(wp, wk):
            m = (a1 + a2) / 2
            stack += [(m, a2, p1, p2, K1, K2), (a1, m, p1, p2, K1, K2)]
        elif wp >= wk:
            m = E.round_down(E.root_bracket(p1 * p2, 2)[0])
            stack += [(a1, a2, m, p2, K1, K2), (a1, a2, p1, m, K1, K2)]
        else:
            m = E.round_down(E.root_bracket(K1 * K2, 2)[0])
            stack += [(a1, a2, p1, p2, m, K2), (a1, a2, p1, p2, K1, m)]
    vol = sum((Fr(b[1]) - Fr(b[0])) * (Fr(b[3]) - Fr(b[2])) * (Fr(b[5]) - Fr(b[4])) for b in boxes)
    require(vol == (A_MAX - A_MIN) * (P_MAX - P_MIN) * (K_MAX - K_MIN), "sub-boxes tile the region")
    res["sub_boxes_certified"] = len(boxes)
    res["sub_box_format"] = ["alpha_lo", "alpha_hi", "p_lo", "p_hi", "K_lo", "K_hi", "j (test point 2^-j)"]
    res["sub_boxes"] = boxes
    res["counters"] = E.counters()
    res["floating_point_used_for_decisions"] = False
    res["status"] = "PASS"
    certlib.write_certificate(out, res)
    print("PASS  C6 existence region: %d sub-boxes" % len(boxes))


if __name__ == "__main__":
    main()
