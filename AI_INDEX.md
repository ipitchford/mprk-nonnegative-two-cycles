# Agent navigation

## Identity and starting point

- Release slug: `mprk-nonnegative-two-cycles`
- Version: `1.0.0-candidate`; DOI 10.5281/zenodo.23511175
- Scholarly creator: Anonymous
- Manuscript: [paper.pdf](paper.pdf) (reference), [paper.md](paper.md) (text), [paper.tex](paper.tex) (source)
- Machine-readable claims and their status: [CLAIMS.json](CLAIMS.json)
- Overview and results table: [README.md](README.md); status: [STATUS.md](STATUS.md)
- Evidence limits: [ASSURANCE.md](ASSURANCE.md); internal review objections without ratings: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)
- Provenance and publication changes: [PROVENANCE.md](PROVENANCE.md)
- Prior-art boundary: [novelty_report.md](novelty_report.md), [SOURCES.md](SOURCES.md), [CITATION_AUDIT.md](CITATION_AUDIT.md); original target: [target.yaml](target.yaml)
- Citation metadata: [CITATION.cff](CITATION.cff), [.zenodo.json](.zenodo.json); licences: [LICENSES.md](LICENSES.md)

## Scientific map

The object is the MPRK22(alpha) scheme applied to positive conservative linear
systems y' = Ay. For two species it reduces to a scalar map F on (0,1)
(paper Section 2, equations (7)-(8)); the steady state is its unique fixed
point p (Lemma 3).

- Theorem 5 (certificates C1, C2): at alpha = 4/7, K = 10^4, p = 1/1000 the
  map has an attracting two-cycle reached from 1/2, so the formalised
  conjecture (Definition 1, Conjecture 2) is false.
- Theorem 6 (C2): complete phase portrait; basin of the steady state is an
  interval of length about 0.108.
- Theorem 9 (proof): two-cycles for every alpha in (1/2, 1/phi) at small p and large K.
- Theorems 10-11 (C6, C4): explicit certified parameter sets.
- Theorem 12 (C3): three-species chain. Theorem 13 (C7): third-order MPRK43II(37/50).
- Proposition 14, Theorems 15-16 (proofs; identities C5): positive results for
  alpha = 1/2 and alpha = 1, two species. Conjecture 17 is open.

Do not infer: failure of MPRK22(1) or MPRK22(1/2) for two species (proved
convergent); necessity of asymmetry or large steps (observed, not proved);
behaviour of nonlinear systems (exploratory only); historical priority beyond
the dated search; formal verification or peer review (neither exists).

## Replay

From the repository root, with Python 3.10 or later and no third-party packages:

```sh
python3 -S -O verify_all.py              # C1-C7, byte-identical (about 15 minutes)
python3 -S -O verify_all.py --skip-slow  # omit C4 and C6
python3 -S -B negative_controls.py       # the checks reject wrong inputs
shasum -a 256 -c SHA256SUMS
shasum -a 256 -c MANIFEST.sha256
```

Entry points: [verify_all.py](verify_all.py), [negative_controls.py](negative_controls.py);
archived certificates in [certificates/](certificates/); exact arithmetic in
[mprk_exact.py](mprk_exact.py) and the maps in [mprk_maps.py](mprk_maps.py).
Exploratory floating-point scripts are in [exploratory/](exploratory/).
Recorded results: [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md); environment: [ENVIRONMENT.txt](ENVIRONMENT.txt).

## Useful next research questions

Prove or refute two-species global convergence of MPRK22(alpha) for
alpha in [1/phi, 1) and alpha > 1 (Conjecture 17(a)); prove global convergence
of MPRK22(1) in every dimension (Conjecture 17(b)); find closed-form thresholds
p0(alpha), K0(alpha); test whether other nonnegative-tableau MPRK and SSP-MPRK
schemes admit attracting cycles on asymmetric stiff systems.
