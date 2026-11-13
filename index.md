# Paper 9 — Solar Capture Tests of Inelastic Dark Matter after LZ

> Reproduction of **Di Mauro & Shaikh (2026), arXiv:2609.06760**
> (`papers/09_Solar_Capture_Inelastic_DM_LZ.pdf`) built with
> **ASTRA** (`astra.yaml`) + **lightcone-cli** (`lc`).
> Every figure and decision below is resolved from `astra.yaml` via MySTRA
> — no values are hard-typed in this document. Preview with `myst start`.

## The question

Can IceCube solar-neutrino limits rule out thermal-Higgsino explanations of
the LZ 248 keV nuclear-recoil candidate?

## Inputs (declared in `astra.yaml`)

- {astra}`input:paper9_pdf` — Paper 9 PDF.
- {astra}`input:lz_eft_paper` — LZ EFT search (Sep 1, 2026): defines the 248 keV event.
- {astra}`input:papers_collection_notes` — `papers/README.md` mapping of all 12 papers.
- {astra}`input:physics_engine` + {astra}`input:bridge_tools` — physics engine and MCP wrappers.
- {astra}`input:paper9_script` — driver producing all four outputs.

## Outputs

:::{astra}output:paper9_solar_capture_vs_splitting
:::

:::{astra}output:paper9_higgsino_exclusion
:::

:::{astra}output:paper9_mixed_dm_constraints
:::

:::{astra}output:paper9_bbn_predictions
:::

All four are produced by one recipe — `{astra}recipe:paper9_solar_capture_vs_splitting`
is representative; see `astra.yaml` per-output `recipe.command`.

Rendered (committed) figures live in `analysis/output/*.png`; materialized
(by-manifest) copies land in `results/<universe>/<output>.<format>` after
`lc materialize local`.

## Decisions

- {astra}`decision:post_capture_cooling`
- {astra}`decision:neff_scenario`
- {astra}`decision:mixed_dm_grid`

Universe in use: `universes/baseline.yaml` — all three cooling scenarios,
full N_eff stack, 60×60 grid (the configuration that produced the committed
`analysis/output/` figures).

## What the figures mean (short)

- **Plot 1 (C_⊙ vs δ):** capture rate falls with endothermic splitting; IceCube
  2Γ_lim bands pick out which δ are still allowed for each model's benchmark.
- **Plot 2 (Γ_A exclusion):** for the Higgsino benchmark δ_LZ ≈ 377 keV lies
  well below δ_lim ≈ 512–565 keV in all three cooling scenarios → thermal
  Higgsino excluded; pseudo-Dirac and O₄ benchmarks are not.
- **Plot 3 (mixed-DM χ² plane):** after Higgsino exclusion the joint χ²
  prefers f_WIMP ~ 0.1–0.3; pure-Higgsino stripe is disfavoured.
- **Plot 4 (BBN vs N_eff):** with a 1.08 TeV WIMP, ΔN_eff,WIMP ≈ 0 while the
  7 keV sterile (sin²2θ ~ 7×10⁻¹¹) gives ΔN_eff ~ 0.01 — BBN markers shift
  only slightly from the SM baseline, comfortably inside Planck's N_eff band.

## Reproducing

```bash
astra validate && astra info
lc status
lc compute launch local   # optional — tune resources
lc materialize local      # each output + manifest committed; byte-identical rebuilds short-circuit
myst start                # preview this page at http://localhost:3000
```

## Toggle

Most documentation for the Hermes fleet is now kept through Lightcone products.
Toggle it off/on without losing anything:

```bash
hermes config set lightcone.enabled false
hermes config set lightcone.enabled true
```
