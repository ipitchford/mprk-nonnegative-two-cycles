# Provenance

## Input

The research package was supplied as an approved Evidence Press candidate
(revised package of 20 October 2026). Its 43 research files matched their
recorded SHA-256 hashes at intake. A private file answering internal reviews
and the private review records were excluded from the public package.

## Research and review

The manuscript, proofs and software were prepared with substantial assistance
from an AI system; the originating model identities are not asserted here.
The original target statement (`target.yaml`) and the dated prior-art search
(`novelty_report.md`) are included.

Internal AI-assisted review comprised a first round by several AI reviewer
roles, a full revision, a verification re-review with independent
recomputation (mpmath interval arithmetic and sympy), and an adversarial review
by the OpenAI model `gpt-5.6-sol` on 20 October 2026. A final round of fixes
was then applied. A rating-free summary of the substantive objections and their
resolution is in [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). This is not external
human peer review.

## Publication changes

For publication the following changes were made to the approved candidate:
the version, DOI and repository were added to the title page, PDF metadata,
reproducibility section and README; one side remark comparing the arXiv and
journal abstracts of the source paper was completed (the journal abstract also
adds three sentences on alpha < 1/2), in both the paper and the novelty report,
after re-checking the Crossref and arXiv abstracts; `negative_controls.py`, the
release records, `paper.md`, the CI workflow and the licence files were added.
The verifiers, libraries and certificates are unchanged; `SHA256SUMS` was
regenerated and `MANIFEST.sha256` covers every shipped file.

Scholarly creator: **Anonymous**. Ian Pitchford is the repository maintainer and
publisher through Evidence Press, not an attributed scholarly author.
