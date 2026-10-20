# Nonnegative Runge–Kutta coefficients do not make modified Patankar schemes globally stable: a counterexample with attracting two-cycles

Research package, 20 October 2026. Author: Anonymous. Unrefereed manuscript with its exact verifiers.

- Version: 1.0.1-candidate (unrefereed research candidate; AI-assisted; not peer reviewed)
- DOI: [10.5281/zenodo.23511175](https://doi.org/10.5281/zenodo.23511175)
- Repository: <https://github.com/ipitchford/mprk-nonnegative-two-cycles>
- Evidence Press page: <https://evidencepress.org/releases/mprk-nonnegative-two-cycles/>
- Licences: prose, certificates and data CC0-1.0 ([LICENSE](LICENSE)); code MIT ([LICENSE-CODE](LICENSE-CODE)); see [LICENSES.md](LICENSES.md)
- Start here if you are an agent: [AI_INDEX.md](AI_INDEX.md); status and evidence limits: [STATUS.md](STATUS.md), [ASSURANCE.md](ASSURANCE.md)

## In plain English

Modified Patankar–Runge–Kutta (MPRK) schemes are time-stepping methods for stiff chemical, biological and environmental models. They never produce negative concentrations and never lose mass, however large the time step. Izgin, Kopecz, Meister and Schilling (Numerical Algorithms 96, 2024; arXiv:2309.01562) conjectured that when the Runge–Kutta coefficients behind such a scheme are all nonnegative, the scheme always ends up at the correct steady state, from any starting point and with any step size.

This package shows that the conjecture is false. Take the scheme MPRK22(4/7), whose nonzero coefficients 4/7, 1/8 and 7/8 are all positive, and a simple two-species linear reaction with rates 9990 and 10. With step size 1, started from (1/2, 1/2), the scheme does not settle at the true steady state (0.001, 0.999). Instead it jumps for ever between two wrong states, about (0.4935, 0.5065) and (1.08e-7, 0.99999989). Positivity and mass are preserved exactly, so the usual sanity checks do not notice. The true steady state is still the only steady state of the scheme and is locally stable, but only starting points in a certified interval of length about 0.108 reach it (about 11% of starting mass fractions, if these are drawn uniformly); we determine exactly which ones.

Every failure we found involves a very lopsided equilibrium and a step much longer than the fastest time scale, and an asymptotic analysis explains why. The original experiments used a symmetric equilibrium, which is consistent with their not seeing it. We show the same behaviour for a whole interval of methods, for a three-species reaction chain and for a third-order scheme. We also prove that two commonly used members, MPRK22(1) and MPRK22(1/2), always converge for two-species systems. The paper ends with advice for users and a refined conjecture.

## What is proved, what is certified, what is exploratory

| Result in the paper | Status | Evidence |
|---|---|---|
| Lemma 3: unique positive fixed point and local stability (two species) | written proof | paper |
| Theorem 5: the counterexample (attracting two-cycle of MPRK22(4/7)) | certified | C1, C2 |
| Theorem 6: complete phase portrait of the counterexample | certified + written proof | C2 |
| Theorem 9: two-cycles for every alpha in (1/2, 1/phi) at small p and large K | written proof | paper |
| Theorem 10: explicit (alpha, p, K) region with two-cycles | certified + short written argument | C6 |
| Theorem 11: attracting two-cycles for every K >= 10^4 at nine (alpha, p) pairs; explicit box around the counterexample | certified | C4 |
| Theorem 12: attracting two-cycle for the chain A <-> B <-> C | certified | C3 |
| Theorem 13: attracting two-cycle for the third-order scheme MPRK43II(37/50) | certified | C7 |
| Proposition 14: no spurious positive fixed points in any dimension, 1/2 <= alpha <= 1 | written proof | paper |
| Theorem 15: global convergence for alpha = 1/2 (two species) | written proof; identities certified | C5 |
| Theorem 16: global convergence for alpha = 1 (two species) | written proof; identities certified | C5 |
| Phase diagram, random networks, error estimator, nonlinear and other-scheme tests | exploratory (floating point) | `exploratory/` |

"Certified" means an exact rational interval-arithmetic computation (Python standard library only) that writes a JSON certificate. "Exploratory" means floating-point numerics without proof.

## Files

- `paper.tex`, `paper.pdf`: the manuscript.
- `mprk_exact.py`: outward-rounded rational interval arithmetic (standard library only).
- `mprk_maps.py`: interval enclosures of the MPRK22 maps; `certlib.py`: input hashes and deterministic JSON output.
- Verifiers: `verify_mprk22_cycle.py` (C1, the original verifier, unchanged), `verify_phase_portrait.py` (C2), `verify_three_state_cycle.py` (C3), `verify_effective_region.py` (C4), `verify_algebraic_identities.py` (C5), `verify_existence_region.py` (C6), `verify_mprk43_cycle.py` (C7).
- `certificates/`: archived JSON certificates C1–C7. Each records SHA-256 hashes of the verifier and library sources, the parameters and all exact rational enclosures.
- `verify_all.py`: replays every certificate in a temporary directory and checks byte identity.
- `negative_controls.py`: checks that the tests reject wrong inputs (symmetric system, alpha = 1/2, 2/3, 1, a tampered certificate or library, the misprint a31 = 113/222) and that no acceptance test uses `assert`.
- `illustrate_mprk22_cycle.py`, `trajectory.csv`, `fig_trajectory.png`, `fig_portrait.png`: illustrations (NumPy, Matplotlib; not proofs).
- `exploratory/`: floating-point scripts X1–X6 and their CSV/JSON/PNG outputs.
- `novelty_report.md`: record of the prior-art search.
- `target.yaml`: the originally frozen target statement.
- `SHA256SUMS`: SHA-256 hashes of the research files listed above.
- Release records: `STATUS.md`, `ASSURANCE.md`, `PROVENANCE.md`, `AI_INDEX.md`, `CLAIMS.json`, `REVIEW_SUMMARY.md`, `RESEARCH_GATES.json` (prior-art, limitation and extension records), `REPLAY_RECEIPT.md`, `SOURCES.md`, `CITATION_AUDIT.md`, `ENVIRONMENT.txt`, `CITATION.cff`, `.zenodo.json`, licence files, `paper.md` (accessible text of the paper) and `MANIFEST.sha256` (hashes of every shipped file).

## How to replay

With Python 3.10 or later and no third-party packages, from this directory:

```bash
python3 -S -O verify_all.py              # all seven certificates (about 15 minutes)
python3 -S -O verify_all.py --skip-slow  # omit C4 and C6
python3 -S -B negative_controls.py       # the checks reject wrong inputs (seconds)
shasum -a 256 -c SHA256SUMS              # file integrity
```

`verify_all.py` copies the verifiers into a fresh temporary directory, runs each with `python -S -O`, and compares each regenerated certificate byte for byte with `certificates/`. A single verifier can be run as

```bash
python3 -S -O verify_phase_portrait.py --out /tmp/new_C2.json
```

The verifiers refuse to overwrite an existing file. Acceptance tests raise an exception instead of using `assert`, so they stay active under `-O`. No floating-point value enters any decision; decimal fields in certificates are labelled as displays only.

Tested with Python 3.13.5 on macOS. The exploratory scripts need NumPy, Matplotlib and mpmath (run them from inside `exploratory/`). The paper builds with `pdflatex paper.tex`, run twice.

## Trust boundary

The verifiers are short Python programs, not proof-assistant developments; correct execution by the Python interpreter is part of the trust boundary. The reduction of each scheme to the evaluated map is written out in the paper. MPRK43II is defined in full in Section 6.3 of the paper, following Kopecz and Meister (BIT 58, 2018; arXiv:1703.05052, eqs. (63) and (68)); as a cross-check the code agrees with PositiveIntegrators.jl (`src/mprk.jl` at commit 37efc009f8baed988b999f56ee0dc76021bf8ee1), and third-order convergence was checked numerically. Hashes and byte-identical replay show integrity (the certificates are what the code produces); whether the code implements the mathematics correctly is checked against the written specifications in the paper and the verifier docstrings. The manuscript was prepared with substantial assistance from an AI system and has not been peer reviewed.

## Main source

Izgin, T., Kopecz, S., Meister, A., & Schilling, A. (2024). On the non-global linear stability and spurious fixed points of MPRK schemes with negative RK parameters. *Numerical Algorithms*, 96(3), 1221–1242. https://doi.org/10.1007/s11075-024-01770-7 (arXiv:2309.01562).
