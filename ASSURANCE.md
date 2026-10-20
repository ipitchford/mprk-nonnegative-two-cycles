# Assurance and trust boundaries

## What supports the result

- **Written proofs.** Lemma 3, Theorem 9, Proposition 14 and Theorems 15–16
  are proved in the paper; Theorems 6 and 10 combine a written argument with
  a certificate.
- **Exact certificates C1–C7.** Seven standalone verifiers (Python standard
  library only) evaluate the stated maps in outward-rounded rational interval
  arithmetic and write JSON certificates that record the SHA-256 hashes of their
  own source and of the library files. No floating-point value enters a
  decision. Acceptance tests raise exceptions rather than using `assert`, so
  they remain active under `python -O`.
- **Byte-identical replay.** `verify_all.py` regenerates all seven certificates
  in a fresh temporary directory and compares them byte for byte with
  `certificates/`. During preparation a separate clean replay of C1–C7 in the
  producer workflow matched byte for byte (C4 about 13 minutes, C6 about 70 s),
  and the publication gate repeats the replay on the frozen archive; see
  [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md).
- **Negative controls.** `negative_controls.py` checks that the trapping test
  holds at the counterexample but fails for the symmetric system and for
  alpha = 1/2, 2/3 and 1; that tampered certificates and libraries are detected;
  that the verifiers stop on a wrong problem or on the misprint a31 = 113/222;
  and that no acceptance test uses `assert`.
- **AI-assisted review.** The manuscript went through internal AI-assisted
  review (several AI reviewer roles, a verification re-review with independent
  recomputations in different software, and an adversarial review by a
  different AI model). The substantive objections and how they were resolved
  are summarised without ratings in [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).
  This is internal scrutiny, not peer review.

## What is not established

- No journal or editorial peer review, no external specialist review and no
  independent reproduction by an unaffiliated party. Recomputations by AI
  reviewers inside the producer workflow are implementation diversity, not
  independent reproduction.
- No proof-assistant formalisation. Correct execution of the verifiers by the
  Python interpreter is part of the trust boundary. Hashes and byte-identical
  regeneration show integrity (the certificates are what the code produces),
  not that the code implements the mathematics; that correspondence is checked
  against the written specifications in the paper and the verifier docstrings.
- The conjecture was stated in words by its authors; the paper refutes its own
  formalisation (Definition 1, Conjecture 2), which it argues is the natural
  literal reading.
- Novelty is bounded by the dated search in [novelty_report.md](novelty_report.md).
- Exploratory results (`exploratory/`) are floating-point experiments, not proofs.
