# Typesetting

- `paper.pdf` is built from `paper.tex` with pdfTeX (TeX Live 2026, pdfTeX 1.40.29),
  running `pdflatex paper.tex` three times in a directory that also contains
  `fig_trajectory.png`, `fig_portrait.png` and `exploratory/X1_phase_diagram.png`.
  Fonts are embedded Latin Modern Type 1 fonts.
- `paper.md` is an accessible text version generated from `paper.tex` with
  Pandoc 3.9 after replacing `\cite`, `\ref` and `\eqref` by the numbers of the
  LaTeX build. The PDF is the reference version.
- The figures are produced by `illustrate_mprk22_cycle.py` and
  `exploratory/X1_phase_diagram.py` (floating point; illustrations, not proofs).
