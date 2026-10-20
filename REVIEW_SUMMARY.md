# Summary of internal review objections and their resolution

This is a public summary of the internal AI-assisted review of this candidate.
It records the substantive scientific objections and what was done about them.
It omits reviewer verdicts and any scores; the original review records are
kept privately. None of this is journal peer review, external specialist
review or independent reproduction by an unaffiliated party.

## What was reviewed and how

1. **First round.** Several AI reviewer roles (editorial, rigorous numerics,
   MPRK schemes, discrete dynamics, and a role arguing against the core claim)
   read an earlier, shorter version. They replayed the verifier, and some
   recomputed the scalar reduction, the cycle and the local multiplier in
   independent high-precision code. No mathematical error was reported.
2. **Revision.** The paper was substantially extended in response (see below).
3. **Verification re-review.** An AI reviewer checked each response, re-derived
   the asymptotics and identities with sympy, recomputed the phase portrait,
   the three-species trapping argument and the MPRK43II scheme with its own
   mpmath interval code (independent of the package library), resolved the
   cited DOIs against Crossref, and replayed all seven certificates.
4. **Adversarial review.** A different AI model (OpenAI `gpt-5.6-sol`,
   20 October 2026) was asked for the strongest counterexample or gap.
5. **Final fixes** were applied and spot-checked, and the full replay of
   C1–C7 was repeated in a clean environment.

## Principal objections and dispositions

| Objection | Disposition |
|---|---|
| The conjecture was paraphrased, not quoted, and the source's notion of "global" was not stated. | The source sentence and its description of the desired property are quoted with page references; Definition 1 and Conjecture 2 are explicitly labelled as the paper's formalisation. |
| The source states the conjecture informally; the paper should not imply the authors stated the formal quantifiers. | Reworded: the paper refutes the natural literal reading, which it states. |
| The paper did not say where the failure occurs or what users should conclude. | Added a mechanism section, an explicit parameter interval (Theorem 9), certified regions (Theorems 10–11), an exploratory phase diagram, and a section for users (adaptive step control, fixed-step risks). |
| Only two-species dynamics were treated. | Added a certified three-species chain (Theorem 12), an exploratory random-network study and a certified third-order example (Theorem 13). |
| Qualitative results lacked explicit constants and attractivity. | Theorems 10 and 11 give explicit certified parameter sets with attracting cycles; closed-form thresholds remain open and are stated as such. |
| Which initial values fail was not determined. | Theorem 6 gives the complete phase portrait of the counterexample with a certified basin interval. |
| Positive results were missing. | Proposition 14 and Theorems 15–16 (global convergence for alpha = 1/2 and alpha = 1, two species) were added, and a refined conjecture stated. |
| Literature was thin. | References were expanded to the MPRK and Patankar literature, stability theory, Izgin's thesis and the classical spurious-dynamics literature; each was checked against Crossref or arXiv. |
| The MPRK43II tableau entry was printed as a31 = 113/222 (found by the re-review and the adversarial review). | Corrected to 73/222 in the theorem and the verifier docstring; the verifier and certificate C7 always used 73/222. C7 was regenerated because it records the verifier hash, and a negative control now checks that the misprint is rejected. |
| MPRK43II was not fully defined in the paper and its external reference implementation was not pinned. | Section 6.3 now defines the scheme in full, citing the source equations, and records the cross-checked implementation commit. |
| Several statements were stronger than the evidence (necessity of asymmetry and large steps, "exact" cycles in floating-point scans, the basin boundary points, the "about 89%" figure). | Reworded as observations of the scans, numerically detected cycles, the repelling boundary cycle, and a measure stated for uniformly drawn mass fractions. |
| The "only if" direction of the refined conjecture for 0 < alpha < 1/2 was attributed to the wrong result. | Attributed to the local multiplier formula of Lemma 3(b), with the range alpha > 0 stated. |
| The comparison of the arXiv and journal abstracts of the source was inaccurate. | Corrected; at publication the remark was completed to note the journal abstract's three further sentences. |
| Signposting for non-specialist readers. | Added a definition of general MPRK schemes, orienting sentences and a clarification that local stability is meant within the invariant mass set. |

## Residual limitations acknowledged

Closed-form thresholds p0(alpha), K0(alpha); a sharp failure region; global
convergence of MPRK22(1/2) for three or more species; part (b) of the refined
conjecture; uniqueness of the positive fixed point of MPRK43II(37/50); and a
replay in an independent interval library such as Arb or MPFI. The paper states
these as open or not claimed.
