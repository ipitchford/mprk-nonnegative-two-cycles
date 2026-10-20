#!/usr/bin/env python3
"""Certificate C3: an attracting two-cycle of MPRK22(4/7) for a three-species chain.

System (paper: theorem 'Three species'): the reversible linear chain A <-> B <-> C with
rate constants
    A -> B : 10,   B -> A : 10^4,   B -> C : 10^4,   C -> B : 1,
i.e. y' = A y with the column-conservative Metzler matrix
    A = [[-10, 10000, 0], [10, -20000, 1], [0, 10000, -1]],
step size h = 1 and alpha = 4/7.  Mass is conserved, so the map is studied in
the coordinates x = (y_A, y_B) with y_C = 1 - y_A - y_B.

The program certifies, in exact outward-rounded rational interval arithmetic,
for two parallelograms B_a = a + S_a [-1,1]^2 and B_b = b + S_b [-1,1]^2:
  (i)   B_a and B_b lie in the open simplex and are disjoint;
  (ii)  G(B_a) is inside the interior of B_b and G(B_b) inside the interior of
        B_a (mean-value enclosures with interval Jacobians);
  (iii) G^2 is a contraction on B_a in the norm ||S_a^{-1} x||_inf, with an
        explicit factor below one;
  (iv)  the equilibrium pi = (1000, 1, 10000)/11001 is a fixed point whose
        2x2 Jacobian satisfies the Jury conditions |det| < 1, |tr| < 1 + det;
  (v)   the exact orbit of y0 = (1/3, 1/3, 1/3) enters B_a or B_b after an
        explicit number of steps (interval enclosure of the orbit).
Uniqueness of the positive fixed point follows from the paper's proposition 'No spurious positive fixed points' (written
proof, 1/2 <= alpha <= 1).  The data a, b, S_a were computed beforehand in
high-precision floating point; they are certificate inputs, and any values
passing the checks would do.
Run:  python -S -O verify_three_state_cycle.py --out NEW.json
"""
from fractions import Fraction as Fr
import mprk_exact as E
from mprk_exact import Iv, require
import mprk_maps as Mm
import certlib

ALPHA = Fr(4, 7)
H = Fr(1)
A = [[-10, 10000, 0], [10, -20000, 1], [0, 10000, -1]]
SYS = Mm.LinearSystem(ALPHA, A, H)

# certificate data (rational; from a prior high-precision computation)
a_pt = [Fr("0.03587363173777860102194947"), Fr("0.3308531780596119557532576")]
b_pt = [Fr("0.1776129707353341553494469"), Fr("7.2405098329874542477599e-10")]
V = [[Fr("0.2982123468478110346"), Fr("-0.1261132220897542396")],
     [Fr("0.9544995527434891617"), Fr("1.087070724072691356")]]
DELTA = Fr(1, 10 ** 5)        # size of B_a
THETA = Fr(3, 2)              # inflation of the linearised image defining B_b


def full(x):
    x = [E.as_iv(v) for v in x]
    return [x[0], x[1], 1 - x[0] - x[1]]


def G2(x):
    y = SYS.step(full(x))
    return [y[0], y[1]]


def DG2(x):
    J = SYS.jacobian(full(x))
    return [[J[i][j] - J[i][2] for j in range(2)] for i in range(2)]


def matmul(X, Y):
    return [[sum((E.as_iv(X[i][k]) * E.as_iv(Y[k][j]) for k in range(2)), Iv.point(0))
             for j in range(2)] for i in range(2)]


def inv2(S):
    (p, q), (r, s) = S
    det = p * s - q * r
    require(det != 0)
    return [[s / det, -q / det], [-r / det, p / det]]


def hull_box(c, S):
    """Axis-aligned interval hull of c + S [-1,1]^2."""
    return [Iv(c[i] - abs(S[i][0]) - abs(S[i][1]), c[i] + abs(S[i][0]) + abs(S[i][1])) for i in range(2)]


def image_in(c_from, S_from, c_to, S_to):
    """Coordinates of G(c_from + S_from [-1,1]^2) relative to the target box.

    Returns S_to^{-1} (G(c) - c_to) + S_to^{-1} DG(hull) S_from [-1,1]^2 as two intervals.
    """
    Sinv = inv2(S_to)
    Gc = G2(c_from)
    J = DG2(hull_box(c_from, S_from))
    off = [Gc[i] - c_to[i] for i in range(2)]
    base = [Sinv[i][0] * off[0] + Sinv[i][1] * off[1] for i in range(2)]
    T = matmul(Sinv, matmul(J, S_from))
    unit = Iv(-1, 1, exact=True)
    return [base[i] + T[i][0] * unit + T[i][1] * unit for i in range(2)], T


def jury(J):
    tr = J[0][0] + J[1][1]
    det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
    return tr, det


def main():
    out = certlib.out_path()
    E.reset_counters()
    params = {"alpha": "4/7", "h": "1", "A": [[str(x) for x in row] for row in A],
              "reaction_chain": "A->B 10, B->A 10000, B->C 10000, C->B 1"}
    res = {"certificate": "C3 attracting two-cycle for a three-species chain",
           "status": "RUNNING", "parameters": params,
           "parameters_sha256": certlib.parameters_hash(params),
           "input_sha256": certlib.input_hashes("verify_three_state_cycle.py")}

    S_a = [[DELTA * V[i][j] for j in range(2)] for i in range(2)]
    Ja = DG2(a_pt)
    S_b_iv = matmul(Ja, S_a)
    # rational S_b: inflate the midpoint of the linearised image
    S_b = [[THETA * S_b_iv[i][j].mid() for j in range(2)] for i in range(2)]
    S_b = [[E.round_up(x) if x > 0 else E.round_down(x) for x in row] for row in S_b]

    # (i) positivity and disjointness
    Ha, Hb = hull_box(a_pt, S_a), hull_box(b_pt, S_b)
    for Hx in (Ha, Hb):
        require(Hx[0].lo > 0 and Hx[1].lo > 0 and (1 - Hx[0] - Hx[1]).lo > 0, "box in open simplex")
    require(Ha[1].lo > Hb[1].hi, "disjoint boxes (second coordinate)")

    # (ii) mutual inclusion
    ab, T_ab = image_in(a_pt, S_a, b_pt, S_b)
    ba, T_ba = image_in(b_pt, S_b, a_pt, S_a)
    for c in ab + ba:
        require(c.lo > -1 and c.hi < 1, "image strictly inside target box")

    # (iii) contraction of G^2 on B_a in the S_a-norm: S_a^{-1} DG(B_b) DG(B_a) S_a = T_ba T_ab
    Mprod = matmul(T_ba, T_ab)
    rows = [max(abs(Mprod[i][0].lo), abs(Mprod[i][0].hi)) + max(abs(Mprod[i][1].lo), abs(Mprod[i][1].hi))
            for i in range(2)]
    contraction = max(rows)
    require(contraction < Fr(3, 4), "G^2 contraction factor below 3/4")

    # (iv) equilibrium and local stability
    pi = [Fr(1000, 11001), Fr(1, 11001), Fr(10000, 11001)]
    for i in range(3):
        require(sum(A[i][j] * pi[j] for j in range(3)) == 0, "A pi = 0")
    Gpi = SYS.step(pi)
    for i in range(3):
        require(Gpi[i].contains(pi[i]), "pi is a fixed point")
    Jp = DG2(pi[:2])
    tr, det = jury(Jp)
    require(abs(det.lo) < 1 and abs(det.hi) < 1, "Jury |det| < 1")
    require(tr.hi < 1 + det.lo and -tr.lo < 1 + det.lo, "Jury |tr| < 1 + det")

    # (v) orbit of the uniform initial state enters a box
    Sa_inv, Sb_inv = inv2(S_a), inv2(S_b)
    def in_box(x, c, Sinv):
        z = [Sinv[i][0] * (x[0] - c[0]) + Sinv[i][1] * (x[1] - c[1]) for i in range(2)]
        return all(zz.lo > -1 and zz.hi < 1 for zz in z)
    y = [Iv.point(Fr(1, 3)), Iv.point(Fr(1, 3)), Iv.point(Fr(1, 3))]
    entered = None
    for n in range(1, 401):
        y = SYS.step(y)
        x = [y[0], y[1]]
        if in_box(x, a_pt, Sa_inv):
            entered = (n, "B_a")
            break
        if in_box(x, b_pt, Sb_inv):
            entered = (n, "B_b")
            break
    require(entered is not None, "orbit of (1/3,1/3,1/3) enters a trapping box")

    rec = lambda X: [[X[i][j].record() for j in range(2)] for i in range(2)]
    res["boxes"] = {
        "B_a": {"centre": [str(v) for v in a_pt], "S_a": [[str(v) for v in row] for row in S_a],
                "hull": [Ha[0].record(), Ha[1].record()]},
        "B_b": {"centre": [str(v) for v in b_pt], "S_b": [[str(v) for v in row] for row in S_b],
                "hull": [Hb[0].record(), Hb[1].record()]}}
    res["inclusions"] = {"G(B_a)_in_B_b_coordinates": [c.record() for c in ab],
                         "G(B_b)_in_B_a_coordinates": [c.record() for c in ba]}
    res["contraction"] = {"norm": "||S_a^{-1} x||_inf", "factor_upper_bound": str(contraction),
                          "factor_display_only": float(contraction),
                          "matrix_enclosure": rec(Mprod)}
    res["equilibrium"] = {"pi": [str(v) for v in pi], "jacobian_trace": tr.record(),
                          "jacobian_det": det.record(),
                          "unique_positive_fixed_point": "proposition No spurious positive fixed points (written proof)"}
    res["initial_state"] = {"y0": ["1/3", "1/3", "1/3"], "enters_box_after_steps": entered[0],
                            "box": entered[1]}
    res["counters"] = E.counters()
    res["floating_point_used_for_decisions"] = False
    res["status"] = "PASS"
    certlib.write_certificate(out, res)
    print("PASS  C3 three-species two-cycle; contraction <= %.6g; y0 enters %s after %d steps"
          % (float(contraction), entered[1], entered[0]))


if __name__ == "__main__":
    main()
