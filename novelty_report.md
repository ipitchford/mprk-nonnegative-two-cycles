# Prior-art search record: positive-coefficient MPRK two-cycles

The original search (including formula-fingerprint searches) was run on 19 October 2026; the status of the conjecture was re-checked and all cited references were verified on 20 October 2026. This record states what was searched and found. It does not establish absolute priority.

## 1. The claim checked

Izgin, Kopecz, Meister and Schilling (Numer. Algorithms 96(3):1221–1242, 2024, https://doi.org/10.1007/s11075-024-01770-7; arXiv:2309.01562v1) conjecture in their abstract, informally, that the linear stability of MPRK schemes "is indeed global if the RK parameters are nonnegative". The desired property is spelled out in their Section 1 (arXiv p. 3): the steady state is a stable fixed point for all step sizes, and the iterates tend to it for all initial values with the same linear invariants. Their evidence for nonnegative tableaux is MPRK22(alpha) on a symmetric two-species system (a = 200, Delta t = 1; Section 3, Figure 3). Their Section 4 (p. 11) notes that for alpha near 1 no restriction on the initial value was found. Izgin's thesis (arXiv:2402.13788v3, Chapter 6) describes investigating or proving the claim as of high importance.

The journal version is open access (CC BY 4.0). Its full text could not be retrieved automatically (publisher login redirect); the abstract was compared through the publisher's Crossref metadata. The conjecture sentence differs from the arXiv version by one comma; elsewhere the two abstracts differ in further punctuation, and the journal abstract adds three sentences on MPRK22(alpha) with alpha < 1/2 (re-checked against the Crossref and arXiv abstracts on 20 October 2026).

## 2. What the package contributes relative to that literature

- A counterexample to the conjecture by an attracting period-two orbit of a standard positive-tableau scheme, with a unique and locally stable physical fixed point (the earlier non-global examples converge to spurious fixed points and have negative tableau entries).
- The complete phase portrait of the counterexample, an asymptotic parameter interval, explicit certified parameter sets, a three-species example and a third-order example (MPRK43II(37/50)).
- Two-species global convergence proofs for alpha = 1/2 and alpha = 1, and a refined conjecture.

## 3. Searches and sources inspected

| Topic | Queries (representative) | Sources inspected | Finding |
|---|---|---|---|
| Original claim and its status | the paper title; "MPRK" "global" "conjecture"; "Patankar" "nonnegative" "global stability" 2025 2026 | arXiv:2309.01562v1 full text; Crossref record of the journal version; Izgin thesis arXiv:2402.13788v3 | Conjecture, definition and experimental design verified; thesis lists the question as open |
| Earlier failures of global stability | "MPDeC(8)" local stability; "MPRK22" negative parameters spurious | Izgin & Öffner, ESAIM M2AN 57(4), 2023; the 2024 paper | Known failures use negative tableau entries and spurious fixed points |
| Periodic or oscillatory behaviour of Patankar schemes | "Patankar" "two-cycle", "period two", "periodic", "oscillation" | Torlo, Öffner & Ranocha, Appl. Numer. Math. 182, 2022; Izgin, Öffner & Torlo, SEMA SIMAI Springer Series, 2024 | Transient oscillations and order reduction studied; no persistent periodic attractor of a positive-tableau scheme found |
| Stability theory of MPRK-type schemes | centre-manifold stability MPRK; SSP-MPRK stability; nonlinear PDS stability | Izgin, Kopecz & Meister (ESAIM M2AN 56, 2022; SIAM J. Numer. Anal. 60, 2022; PAMM 22, 2023); Huang, Izgin, Kopecz, Meister & Shu (ESAIM M2AN 57, 2023) | All local results; no global result for positive tableaux found |
| Step-size control and stability | step size controllers MPRK stability 2026 | Izgin & Ranocha, J. Sci. Comput. 109(3), 2026 | Practical stability concern; not a resolution of the conjecture |
| Classical spurious dynamics | spurious solutions, period-two orbits, Runge–Kutta equilibria | Iserles 1990; Hairer, Iserles & Sanz-Serna 1990; Iserles, Peplow & Stuart 1991; Yee, Sweby & Griffiths 1991; Humphries 1993; Stuart & Humphries 1996 | General phenomenon is classical for nonlinear ODEs; the present mechanism (linear ODE, positive tableau, Patankar weights) is not covered there |
| Formula fingerprints (19 October) | the scalar map, the limit map, the exponent condition r(r-1) > 1, the golden-ratio threshold | web and literature search after each reduction | No match found |
| Status re-check (20 October) | MPRK global stability conjecture nonnegative parameters counterexample periodic | web search; arXiv listing; Kassel institutional repository record of the 2024 paper | No resolution of the conjecture found |

All bibliographic details cited in the paper were checked against Crossref records resolved from their DOIs, or against arXiv metadata, on 20 October 2026. The Stuart–Humphries monograph has no DOI; its details were checked against a Crossref-indexed book review.

## 4. Limits of the search

No author enquiry was made. No paid citation index was used. Non-English and dissertation searches were not exhaustive. A negative web search cannot exclude unpublished or very recent work. Supporting results (Lemma 3, Proposition 14, the alpha = 1 bound) may have earlier proofs in other formulations; the paper claims them as self-contained supporting results, not as priority claims.
