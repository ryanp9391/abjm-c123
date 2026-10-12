# abjm-c123

Research repo (Paul Ryan) computing tree-level structure constants in ABJM theory, in the full
su(4) sector.

**Check `NEXT.md` first** — where the last session left off.

## Layout

- `Paul/Mathematica/su4-alternating/` — the su(4) alternating spin chain (fundamental ⊗ anti-fundamental ⊗ …):
  SoV basis in the companion twist frame, SoV localisation, rotations acting on few sites with free parameters
  (à la arXiv:2603.26617). **Start at its `README.md`** (conventions, naming, results, open items, file map).
  Moved here from the gaudin-cg repo on 2026-10-11.
- `Paul/Mathematica/SU3_Tests/` — earlier su(3) tests (`L2.wb`, `L3_partial.wb`).
- `logbook.tex` — research logbook, maintained by the `logbook-agent` ("update the logbook").

## Conventions (su4-alternating thread)

- Lax ℒ(u) = Σ e_ij ⊗ ((u−θ)δ_ij + h ρ(E_ij)), monodromy with site 1 leftmost, twist on the right, companion
  twist G_c(z). In this Lax rep λ is Ryan–Volin's λ*.
- **Naming (Paul):** "fundamental" site = SoV layer built with t₁ⁿ (ρ(E_ij) = −e_ji); "anti-fundamental" site =
  built with t_n (ρ(E_ij) = e_ij). The L = 2 notebook uses the opposite names — see the README.
- Finite generic twist; no untwisting / small-ε expansions unless Paul asks.

## Notebooks

- **All computations and checks go in wolfbook `.wb` notebooks — never python.** Python only if the wolfbook kernel is
  blocked (say so). The `checks/*.py` scripts are legacy; port to a notebook when needed instead of extending them.
- `.wb` notebooks are wolfbook notebooks; use the wolfbook MCP tools, not generic JSON edits, when VS Code has them
  open. Notebooks built by agents as JSON have **no native outputs** (results in markdown cells) — rerun in VS Code.
- Multi-statement Wolfram cells need `;` between statements; a comment containing `m*)` closes early; memoised
  definitions cache bad results — rerun downstream after a fix.
