# Status

Version: 1.0.0-candidate. DOI: 10.5281/zenodo.23511175. Scholarly creator: Anonymous.

This is an AI-assisted, unrefereed research candidate. It has not been peer
reviewed by a journal, and no external specialist review or independent
reproduction by an unaffiliated party is claimed.

Principal claim: the natural universal reading of the conjecture of Izgin,
Kopecz, Meister and Schilling (Numer. Algorithms 96, 2024) that MPRK schemes
with nonnegative Runge–Kutta coefficients are globally linearly stable is
false. MPRK22(4/7), all of whose nonzero tableau coefficients are positive,
has a certified attracting two-cycle on the two-species linear system with
transfer rates 9990 and 10 at step size 1, while the physical steady state
(0.001, 0.999) is the unique positive fixed point and locally asymptotically
stable. The conjecture was stated in words; the formal statement refuted
(Conjecture 2 of the paper) is the paper's own formalisation of it.

Status of the individual results (proof / certified / exploratory) is listed
in [README.md](README.md) and [CLAIMS.json](CLAIMS.json). Certified results
rest on exact rational interval-arithmetic verifiers executed by Python; they
are not proof-assistant formalisations. Exploratory results use floating
point and are not proofs. Novelty is bounded by the dated search in
[novelty_report.md](novelty_report.md), not proved.
