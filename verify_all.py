#!/usr/bin/env python3
"""Replay every exact certificate in a fresh temporary directory and compare the
regenerated JSON byte for byte with the archived copy in certificates/.

Usage:  python3 -S -O verify_all.py            (all seven certificates)
        python3 -S -O verify_all.py --skip-slow (omit C4 and C6, which take minutes)

Each verifier is run as a separate process with `python -S -O`, from inside the
temporary directory, so archived files can never be overwritten.  Standard
library only.  Exit status 0 means every replayed certificate matched.
"""
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIBS = ["mprk_exact.py", "mprk_maps.py", "certlib.py"]
# (label, verifier, archived certificate, slow?)
JOBS = [
    ("C1", "verify_mprk22_cycle.py", "C1_mprk22_cycle_certificate.json", False),
    ("C2", "verify_phase_portrait.py", "C2_phase_portrait.json", False),
    ("C3", "verify_three_state_cycle.py", "C3_three_state_cycle.json", False),
    ("C4", "verify_effective_region.py", "C4_effective_region.json", True),
    ("C5", "verify_algebraic_identities.py", "C5_algebraic_identities.json", False),
    ("C6", "verify_existence_region.py", "C6_existence_region.json", True),
    ("C7", "verify_mprk43_cycle.py", "C7_mprk43_cycle.json", False),
]


def main():
    skip_slow = "--skip-slow" in sys.argv[1:]
    results = []
    with tempfile.TemporaryDirectory(prefix="mprk22-replay-") as tmp:
        tmp = Path(tmp)
        for f in LIBS + [j[1] for j in JOBS]:
            shutil.copy2(HERE / f, tmp / f)
        for label, script, cert, slow in JOBS:
            if slow and skip_slow:
                results.append((label, "SKIPPED", 0.0))
                continue
            t0 = time.time()
            if label == "C1":
                # the original verifier writes its receipt next to itself
                cmd = [sys.executable, "-S", "-O", script]
                produced = tmp / "mprk22_cycle_certificate.json"
            else:
                produced = tmp / ("replay_" + cert)
                cmd = [sys.executable, "-S", "-O", script, "--out", str(produced)]
            proc = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True)
            dt = time.time() - t0
            if proc.returncode != 0:
                results.append((label, "FAIL (verifier error: %s)" % proc.stderr.strip().splitlines()[-1:], dt))
                continue
            same = produced.read_bytes() == (HERE / "certificates" / cert).read_bytes()
            results.append((label, "PASS (byte-identical)" if same else "FAIL (certificate differs)", dt))
            print("%s  %-32s %s  [%.1f s]" % (label, script, results[-1][1], dt), flush=True)
    ok = all(r[1].startswith("PASS") or r[1] == "SKIPPED" for r in results)
    print("\nSummary:")
    for label, status, dt in results:
        print("  %s: %s" % (label, status))
    print("ALL REPLAYED CERTIFICATES MATCH" if ok else "SOME CERTIFICATES FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
