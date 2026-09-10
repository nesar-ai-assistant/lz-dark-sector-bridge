# Papers Collection — LZ Dark Sector Bridge

Curated collection of papers at the intersection of the LZ direct detection result
(Sep 2026) and cosmology: relic abundance, BBN, indirect detection, and dark sector models.

## The LZ Result

| # | File | arXiv | Title | Relevance |
|---|------|-------|-------|-----------|
| 1 | `01_LZ_2026_EFT_Nuclear_Recoil_Search.pdf` | LZ preprint (Sep 1, 2026) | Search for dark matter particle interactions in an extended nuclear recoil energy window with the LUX-ZEPLIN (LZ) experiment | **The paper.** 248 keV event, 2.6σ global, EFT operators ℒ₁–ℒ₂₀, inelastic models. |

## Cosmological Interpretations & Relic Abundance

| # | File | arXiv | Title | Relevance |
|---|------|-------|-------|-----------|
| 2 | `02_DM_Kinematic_Edge_Relic_Abundance.pdf` | 2609.02608 | Dark Matter at the Kinematic Edge: Interpreting the 248 keV LZ Nuclear-Recoil Candidate | Relic abundance + indirect detection constraints on LZ event. Thermal pseudo-Dirac fermion with σ_N ≈ 6.5×10⁻⁴³ cm², δ ≈ 297 keV. Higgsino at 1.1 TeV with δ ≈ 377 keV. |
| 3 | `03_Higgsino_DM_LZ_248keV.pdf` | 2609.01583 | Higgsino Dark Matter Interpretation of the LUX-ZEPLIN 248 keV Nuclear-Recoil Event | Pure Higgsino ~1 TeV: thermal relic abundance fixes mass, off-diagonal Z coupling gives inelastic scattering. δ ~ 350 keV. |
| 4 | `04_TeV_Higgsino_LZ_FermiLAT.pdf` | 2609.01590 | TeV Higgsino Dark Matter from LZ Nuclear Recoil to Fermi-LAT Gamma Rays | Joint LZ + Fermi-LAT interpretation. Freeze-out fixes mass at 1.1 TeV. Sommerfeld-enhanced annihilation. Predicts gamma-ray line at CTAO sensitivity. |
| 5 | `05_PQ_Origin_Inelastic_EW_DM.pdf` | 2609.02807 | A Peccei–Quinn Origin for Inelastic Electroweak Dark Matter after LUX-ZEPLIN | EW doublet as DM subcomponent (ξ_χ ≈ 0.12–0.15). Mixed DM composition — directly relevant to our mixed CDM+WDM framework. |

## Inelastic DM Models (Exothermic/Endothermic)

| # | File | arXiv | Title | Relevance |
|---|------|-------|-------|-----------|
| 6 | `06_Exothermic_DM_at_LZ.pdf` | 2609.05204 | Exothermic Dark Matter at LZ | Down-scattering with δ ~ 0.5–1 MeV. Both freeze-out and freeze-in cosmologies. Predictions for Ar and Ge targets. |
| 7 | `07_Exo_Endo_Inelastic_DM_Sideband.pdf` | 2609.04673 | Exothermic and Endothermic Inelastic Dark Matter Interpretations at LZ: Sideband Constraints and Future Prospects | LZ sideband analysis constrains exothermic models. Predicts ~4.5 events in 1000 live-day exposure. XENONnT/PandaX cross-checks. |

## Multi-Messenger & Complementary Probes

| # | File | arXiv | Title | Relevance |
|---|------|-------|-------|-----------|
| 8 | `08_Multi_Messenger_Paleo_Detector_LZ.pdf` | 2609.06640 | Multi-Messenger and Paleo-Detector Probes of the LZ Dark Matter Signal | Leptophobic inelastic DM with vector mediator. Multi-messenger search strategy. |
| 9 | `09_Solar_Capture_Inelastic_DM_LZ.pdf` | 2609.06760 | Solar Capture Tests of Inelastic Dark Matter after the LZ High-Recoil Event | IceCube solar neutrino constraints **exclude** thermal Higgsino interpretation (δ ≳ 0.51 MeV required, but LZ needs δ ~ 0.38 MeV). Critical paper. |
| 10 | `10_LZ_to_Colliders_Higgsino.pdf` | 2609.08712 | From LUX-ZEPLIN to Colliders: Probing Higgsino Dark Matter | 1.1 TeV Higgsino with δ ≈ 350 keV → chargino mass splitting O(350 MeV). Collider phenomenology. |
| 11 | `11_ALP_Mediated_Inelastic_DM_LZ.pdf` | 2609.08893 | ALP-mediated inelastic dark matter and the LUX-ZEPLIN high-recoil candidate event LZ230616 | Axion-like particle mediator coupled to gluons. Momentum-dependent interaction reshapes spectrum. |

## BBN + WIMP Cosmology

| # | File | arXiv | Title | Relevance |
|---|------|-------|-------|-----------|
| 12 | `12_BBN_WIMP_Freeze_Out_Modified_Gravity.pdf` | 2608.09390 | Big-Bang Nucleosynthesis and WIMP Dark Matter Freeze-Out in Modified Gravity | BBN + freeze-out in Yukawa cosmology. Complementary probes of modified gravity via BBN and relic abundance. Directly relevant to our joint χ² framework. |

---

## Key Themes for the Bridge Server

1. **Relic abundance as discriminator**: Papers 2–5 all use Ω_χ h² to constrain the WIMP mass (fixing it near 1 TeV for thermal Higgsino). Our `predict_wimp_signal` tool does exactly this.

2. **Mixed DM composition**: Paper 5 explicitly considers subcomponent DM (ξ_χ ≈ 0.12–0.15), directly mapping to our `predict_mixed_dm` and `mixed_dm_budget` framework.

3. **BBN as constraint**: Paper 12 shows BBN + freeze-out as complementary probes. Our `chi2_bbn` and `joint_chi2` implement this.

4. **N_eff bridge**: The sterile neutrino ΔN_eff from our sterile-nu-bbn-server adds to any WIMP contribution — the joint N_eff is the key observable connecting direct detection to primordial nucleosynthesis.

5. **Solar capture exclusion**: Paper 9 is critical — it excludes the canonical thermal Higgsino interpretation via IceCube. This motivates considering alternative DM compositions (mixed CDM+WDM).

---

*Last updated: September 9, 2026*
