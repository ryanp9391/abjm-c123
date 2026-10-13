# Writing guide for the complete notes (2026-10-13)

Master file: `abjm-complete-notes.tex`. It `\input`s `macros.tex` and the section files in `sections/`. Each agent writes ONLY its assigned files.

## Audience and style
- **Readers:** Paul Ryan (the research lead) and expert collaborators in integrability and N = 6 / ABJM. They know Bethe ansatz, SoV, Q-systems and Ryan–Volin. Write like a research-notes appendix: precise, dense, no tutorial material.
- **Model:** match the style of the existing `../abjm-transfer-matrices-symmetries.tex`. Read it first; its sections are the backbone and must be carried over, updated where superseded.
- **Statements:** every statement gets an evidence tag and, where it exists, the notebook it lives in, e.g. `[\nb{pair-group-rotation.wb}]`.
- **Formulas:** state results as displayed equations with all constants. Prefer closed forms; give numbers at the standard point only when no closed form exists.
- **Corrections:** do NOT invent results. Use only what is in the source files. If two sources disagree, use the later one and add a footnote naming the earlier statement. For the global list of known corrections see "Superseded" below.
- **Length:** about 3–6 pages per section file.

## LaTeX rules (strict; pdflatex, MiKTeX)
- **Characters:** pure ASCII source. NO Unicode characters at all (no θ, ⟨, ⊗, —, é in source; write `\theta`, `\langle`, `\otimes`, `--`, `\'e`).
- **Macros:** use ONLY `macros.tex` macros, amsmath/amssymb, booktabs, longtable and array. If you need a convenience macro, define it at the top of YOUR file with `\providecommand`, prefixed with your section name (e.g. `\providecommand{\sovW}{...}`), and list it in your report.
- **Evidence tags:** `\ex` = exact (rational/symbolic) in a notebook; `\hp` = high-precision numerics in a notebook; `\num{10^{-12}}` = numerics with that residual; `\npy` = python only.
- **Labels (fixed; use them for cross-references):**
  - Sections: `sec:overview`, `sec:conv`, `sec:tm`, `sec:abjm`, `sec:loc`, `sec:twvac`, `sec:states`, `sec:sym`, `sec:C`, `sec:configs`, `sec:fact`, `sec:polforms`, `sec:sov`, `sec:3pt`, `sec:odd`, `sec:kappa`, `sec:admiss`, `sec:qform`, `sec:open`, `sec:index`, `sec:superseded`.
  - Subsections you create: label them with your section prefix (e.g. `sec:sov-dict`).
- **File format:** each file starts with `\section{...}\label{...}` (several `\section`s per file are fine if assigned). No `\documentclass`, no `\begin{document}`.
- **Compile test:** create `test-<yourfile>.tex` in `complete-notes/` containing
  ```latex
  \documentclass[11pt]{article}
  \input{macros}
  \begin{document}
  \input{sections/<yourfile>}
  \end{document}
  ```
  Run it from `complete-notes/` as `pdflatex -interaction=nonstopmode -halt-on-error test-<yourfile>.tex`, twice. It must compile without errors; undefined references to OTHER sections' labels are fine. Then delete your `test-<yourfile>.*` files.

## Conventions (must match; details in `../../README.md` and `../w-a-reference-status.md`)
- **Chain:** A F A F ..., site 1 leftmost. A sites carry the 4: rho_A(E_ij) = e_ij, the group acts by g. F sites carry the 4bar: rho_F(E_ij) = -e_ji, the group acts by g^{-T}. Bras are row vectors and operators act from the right.
- **Monodromy and transfer matrices:** Lax L_s(u) = sum e_ij (x) ((u - theta_s) delta_ij + h rho_s(E_ij)); T = L_1 ... L_{2L}; twist on the right; t_a(u) = sum of quantum minors (rows antisymmetrised, column k at u - (k-1)h).
- **Twists:** diagonal twist Lambda = diag(z) (physical frame); companion twist G_c(z). Frame map t^{G_c} R = R t^Lambda, with R = V^{-T} (A) (x) V (F), where V has columns (z_j^3, z_j^2, z_j, 1).
- **Vacua and reference bras:** SoV vacuum <0| = (e_1 (x) e_4)^L in the companion frame; vac = e_1 (x) e_4 per pair in the diagonal frame; W_a = (e_1 + a e_4) (x) (e_4 - a e_1).
- **ABJM point:** theta_A - theta_F = 2h, h = 1/2; u_0 = theta_F + h = theta_A - h.

## Superseded / corrected (respect these everywhere)
1. **Even bridge.** "The even bridge is not reachable / the F-site gap" (old note section 6) is SOLVED: the pair kappa-family K_kappa exists at the ABJM point with z_2 = z_3 (sec:kappa).
2. **theta_A = theta_F.** All computations there are superseded (singular point).
3. **SoV at the ABJM point.** "The Ryan–Volin B is degenerate at the ABJM point" is resolved by the regularised basis (d -> 2h limit) of `sov-abjm-point-limit.wb`.
4. **Q-form of <W_a^L|Psi>.** "<W_a^L|Psi> has no Q-form" (morning of 2026-10-13) is WITHDRAWN. Its overlaps equal those of the twisted vacuum with w = (w_1, z_2, z_3, w_4) (sec:qform).
5. **K_kappa on <0|_comp.** "K_kappa fails on <0|_comp" is refined: the five-word K_kappa is specific to the W class. Block-BPS references admit class-dependent words but have zero cross pairings (sec:admiss).
6. **K_1.** "K_{kappa=1}(tau) = 0 for all non-vacuum states" is FALSE: K_1 vanishes only for n <= 1.
7. **Even/odd ratio.** There is NO even/odd ratio relation.
8. **Untwisting of a rotated pair.** It uses D = diag(c - phi_j), not c + phi_j.
9. **Rerun of 2026-10-13** (`../rerun-2026-10-13.md`):
   - `pair-group-rotation.wb` cell 17: four 2-factor pair words localise with constants (51/31 or 71/31) x prod lambda_k, not prod lambda_k. K_kappa is unaffected.
   - `even-bridge-lpairs.wb`: in the section 2 untwisting table, the F1 and F2 labels are swapped relative to the text. The limits are the same.
10. **The SoV norm B.** It is known to Paul; agent E's derivation in Q's is a cross-check only.
