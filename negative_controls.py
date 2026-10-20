#!/usr/bin/env python3
"""Negative controls for the exact verifiers (standard library only).

A replay that cannot fail is not evidence.  This script checks that the
package's checks *reject* plausible wrong inputs, and accept the right one:

  S1-S5  Semantic discrimination.  The trapping test of Theorem 5 (the two
         rational intervals L = [1e-7, 117e-9] and H = [12/25, 507/1000] map
         strictly into one another) is evaluated with the library map of
         mprk_maps.py.  It must hold for the counterexample
         (alpha = 4/7, p = 1/1000, K = 10^4) and must FAIL for the symmetric
         system p = 1/2, for alpha = 1 and alpha = 1/2 (where the paper proves
         global convergence for two species), and for alpha = 2/3 > 1/phi.
  T1     Certificate tampering: a regenerated C2 certificate equals the archived
         one, and a copy with one altered digit is detected.
  T2     Library tampering: appending one comment line to mprk_maps.py changes
         the regenerated C2 certificate (input hashes are bound).
  T3     Wrong problem: verify_phase_portrait.py with p = 1/2 instead of 1/1000
         must stop with a CertificateError.
  T4     Transcription error: verify_mprk43_cycle.py with the misprint
         a31 = 113/222 (instead of 73/222) must stop with a CertificateError.
  A1     No acceptance test relies on `assert` (which `python -O` removes).

Run from the package root:  python3 -S -B negative_controls.py
                       and  python3 -S -B -O negative_controls.py
Exit status 0 means every control behaved as expected.  Nothing in the
package directory is modified; all copies live in a temporary directory.
"""
import ast
import shutil
import subprocess
import sys
import tempfile
from fractions import Fraction as Fr
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import mprk_exact as E  # noqa: E402
import mprk_maps as Mm  # noqa: E402

LIBS = ["mprk_exact.py", "mprk_maps.py", "certlib.py"]
RESULTS = []


def record(name, expected, observed, detail=""):
    ok = expected == observed
    RESULTS.append(ok)
    print("%-3s %-58s expected=%-6s observed=%-6s %s %s" % (
        name.split(":")[0], name.split(":", 1)[1].strip(), expected, observed,
        "OK" if ok else "UNEXPECTED", detail), flush=True)


def trap_holds(alpha, p, K):
    """True iff F(L) lies in the interior of H and F(H) in the interior of L."""
    M = Mm.Method(Fr(alpha))
    kap = Fr(1, 1) / Fr(K)
    low = E.Iv(Fr(1, 10 ** 7), Fr(117, 10 ** 9), exact=True)
    high = E.Iv(Fr(12, 25), Fr(507, 1000), exact=True)
    f_low = Mm.F(M, low, Fr(p), kap)
    f_high = Mm.F(M, high, Fr(p), kap)
    return (high.lo < f_low.lo and f_low.hi < high.hi and
            low.lo < f_high.lo and f_high.hi < low.hi)


def run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def fresh_copy(tmp, name):
    d = Path(tmp) / name
    d.mkdir()
    for f in LIBS + ["verify_phase_portrait.py", "verify_mprk43_cycle.py"]:
        shutil.copy2(HERE / f, d / f)
    return d


def py(*args):
    # Children inherit this interpreter's optimisation flag, so the controls
    # are exercised in the same mode as the caller.
    flags = ["-S", "-B"] + (["-O"] if sys.flags.optimize else [])
    return [sys.executable] + flags + list(args)


def main():
    print("Python optimisation level:", sys.flags.optimize)
    # ---- semantic discrimination -------------------------------------
    record("S1: trap holds at the counterexample (4/7, 1/1000, 1e4)", True,
           trap_holds(Fr(4, 7), Fr(1, 1000), 10 ** 4))
    record("S2: trap fails for the symmetric system p = 1/2", False,
           trap_holds(Fr(4, 7), Fr(1, 2), 10 ** 4))
    record("S3: trap fails for alpha = 1 (Theorem 16)", False,
           trap_holds(Fr(1), Fr(1, 1000), 10 ** 4))
    record("S4: trap fails for alpha = 1/2 (Theorem 15)", False,
           trap_holds(Fr(1, 2), Fr(1, 1000), 10 ** 4))
    record("S5: trap fails for alpha = 2/3 > 1/phi", False,
           trap_holds(Fr(2, 3), Fr(1, 1000), 10 ** 4))

    archived = (HERE / "certificates" / "C2_phase_portrait.json").read_bytes()
    with tempfile.TemporaryDirectory(prefix="mprk22-controls-") as tmp:
        # ---- T1 certificate tampering --------------------------------
        d = fresh_copy(tmp, "clean")
        r = run(py("verify_phase_portrait.py", "--out", "C2.json"), d)
        regenerated = (d / "C2.json").read_bytes() if r.returncode == 0 else b""
        record("T1: clean C2 regenerates byte-identically", True,
               regenerated == archived)
        text = archived.decode()
        i = next(k for k, ch in enumerate(text) if ch.isdigit() and ch != "7")
        tampered = (text[:i] + "7" + text[i + 1:]).encode()
        record("T1: one-digit change to archived C2 is detected", False,
               regenerated == tampered)
        # ---- T2 library tampering ------------------------------------
        d = fresh_copy(tmp, "lib")
        with open(d / "mprk_maps.py", "a") as fh:
            fh.write("# harmless-looking edit\n")
        r = run(py("verify_phase_portrait.py", "--out", "C2.json"), d)
        out = (d / "C2.json").read_bytes() if r.returncode == 0 else b""
        record("T2: edited library changes the regenerated C2", False,
               out == archived)
        # ---- T3 wrong problem ----------------------------------------
        d = fresh_copy(tmp, "sym")
        src = (d / "verify_phase_portrait.py").read_text()
        if src.count("P = Fr(1, 1000)") != 1:
            raise SystemExit("control setup failed: P definition not found once")
        (d / "verify_phase_portrait.py").write_text(src.replace("P = Fr(1, 1000)", "P = Fr(1, 2)"))
        r = run(py("verify_phase_portrait.py", "--out", "C2.json"), d)
        record("T3: C2 verifier accepts p = 1/2", False,
               r.returncode == 0 or "CertificateError" not in r.stderr,
               "(exit %d: %s)" % (r.returncode, r.stderr.strip().splitlines()[-1:]))
        # ---- T4 transcription error ----------------------------------
        d = fresh_copy(tmp, "a31")
        src = (d / "verify_mprk43_cycle.py").read_text()
        old = "A31 = A21 - 1 / (4 * GAMMA)"
        if src.count(old) != 1:
            raise SystemExit("control setup failed: A31 definition not found once")
        (d / "verify_mprk43_cycle.py").write_text(src.replace(old, "A31 = Fr(113, 222)"))
        r = run(py("verify_mprk43_cycle.py", "--out", "C7.json"), d)
        record("T4: C7 verifier accepts a31 = 113/222", False,
               r.returncode == 0 or "CertificateError" not in r.stderr,
               "(exit %d: %s)" % (r.returncode, r.stderr.strip().splitlines()[-1:]))

    # ---- A1 no assert-based acceptance tests ---------------------------
    files = sorted(HERE.glob("verify_*.py")) + [HERE / f for f in LIBS]
    asserts = []
    for f in files:
        for node in ast.walk(ast.parse(f.read_text(), filename=str(f))):
            if isinstance(node, ast.Assert):
                asserts.append("%s:%d" % (f.name, node.lineno))
    record("A1: assert statements in verifiers/libraries", False, bool(asserts),
           ", ".join(asserts))

    ok = all(RESULTS)
    print("ALL NEGATIVE CONTROLS BEHAVED AS EXPECTED" if ok else "SOME CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
