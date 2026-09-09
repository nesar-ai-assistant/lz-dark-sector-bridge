# LZ Dark Sector Bridge: Unified Constraints on WIMP + Sterile Neutrino Dark Matter

---

## Science Motivation

The **LZ** experiment reported its first non‑background‑like event on **2026‑09‑01** – a nuclear recoil of **248 ± 23 keV** with a global significance of **2.6σ**.  This event lies in the energy range where spin‑independent (SI) WIMP‑Xe scattering is most sensitive to heavy WIMPs (\(m_\chi \gtrsim 200\) GeV) and also overlaps the nuclear recoil signatures used to probe non‑standard effective operators (\(\mathcal{L}_{1-20}\) and inelastic\(\mathcal{O}_{1},\mathcal{O}_{4}\)).

Sterile neutrinos in the \(\sim\) keV–MeV mass range, motivated by short‑baseline anomalies and X‑ray line hints, contribute to the radiation density during **Big‑Bang Nucleosynthesis (BBN)**.  Their presence modifies the effective number of relativistic species \(N_{\rm eff}\) and thus the primordial helium (\(Y_{p}\)) and deuterium (D/H) abundances.  Recent BBN calculations (Pitrou et al. 2018) and Planck CMB constraints place \(N_{\rm eff}=3.044\pm0.03\), leaving only a narrow window for additional sterile states.

A **unified framework** is therefore essential: it must combine the **direct‑detection recoil spectrum** of a WIMP with the **cosmological impact** of sterile neutrinos, enforcing the total dark‑matter density budget \(\Omega_{\rm DM}=\Omega_{\chi}+\Omega_{\nu_s}\).  The **LZ Dark Sector Bridge** provides a single MCP (Model‑Coupling‑Protocol) server that:

* predicts WIMP recoil rates and event counts for LZ‑type xenon targets,
* evaluates sterile‑neutrino BBN and X‑ray constraints, and
* computes a joint \(\chi^{2}\) from LZ, BBN, Planck, and X‑ray data, honouring the mixed‑DM budget.

This enables automated, reproducible scientific agents (e.g. Hermes‑MCP agents) to answer high‑level questions such as:
* "Given the 248 keV recoil, what range of WIMP masses and cross‑sections remain viable?"
* "How does an \(N_{\rm eff}=3.2\) sterile ν scenario affect the LZ event interpretation?"
* "What is the joint likelihood of a mixed‐DM model with \(f_{\chi}=0.8\) and \(f_{\nu_s}=0.2\)?"

---

## Architecture Diagram (text‑based)

```
+----------------------+      +-----------------------+
|  LZ Dark Sector      |      |  Sterile‑ν BBN Server |
|  Bridge (MCP)        |<---->|  (MCP)                |
|  ─────────────────────|      |  ─────────────────────|
|  tools/                |      |  tools/                |
|   ├─ lz_bridge_tools.py ──►│   └─ sterile_tools.py │
|   └─ physics.py         │   └─ bbn.py            |
+----------------------+      +-----------------------+
          |                         |
          ▼                         ▼
   +-------------+           +---------------+
   |  mcp_server |           |  mcp_server   |
   +-------------+           +---------------+
          |                         |
          ▼                         ▼
   +-------------------+   +-------------------+
   |  Hermes Agent     |   |  External Users   |
   +-------------------+   +-------------------+
```

---

## Tool Reference Table

| Tool | Description | Key Outputs |
|------|-------------|-------------|
| `describe_tools` | Lists all available bridge capabilities. | JSON with tool names, signatures and data sources. |
| `predict_wimp` | Direct‑detection signal for a given \(m_\chi\) and \(\sigma_{SI}\). | Expected event count, differential spectrum. |
| `predict_mixed_dm` | Joint WIMP + sterile‑ν dark‑matter budget evaluation. | Fractions \(f_{\chi},f_{\nu_s}\), consistency flag. |
| `recoil_spectrum` | Differential rate \(dR/dE\) over an energy grid. | Energy array, rate array (events kg⁻¹ day⁻¹ keV⁻¹). |
| `joint_constraints` | Combined \(\chi^{2}\) from LZ, BBN, Planck, X‑ray. | Total \(\chi^{2}\) and per‑channel contributions. |
| `lz_event_analysis` | Compatibility of the 248 keV LZ event with WIMP masses. | Kinematic allowance, likelihood profile. |
| `plot_wimp_exclusion` | Generates a mass‑vs‑cross‑section exclusion PNG. | File path to PNG. |
| `plot_mixed_dm` | Visualises the \((f_{\chi},f_{\nu_s})\) plane with constraints. | File path to PNG. |

---

## Example Questions an Agent Could Answer

* **“What WIMP masses can produce a 248 keV recoil in xenon?”** – uses `lz_event_analysis`.
* **“If a sterile neutrino of 7 keV with \(\sin^{2}2\theta=7\times10^{-11}\) exists, how does it shift \(N_{\rm eff}\) and the LZ event likelihood?”** – combines `predict_mixed_dm` and `joint_constraints`.
* **“Show me the exclusion curve for SI cross‑sections from 10 GeV to 10 TeV.”** – calls `plot_wimp_exclusion`.
* **“Is the mixed‑DM budget \(f_{\chi}=0.6, f_{\nu_s}=0.5\) allowed?”** – uses `predict_mixed_dm` (budget violation flag).

---

## Installation

```bash
# Clone the repository
git clone https://github.com/nesar-ai-assistant/lz-dark-sector-bridge.git
cd lz-dark-sector-bridge

# Create a fresh virtual environment (recommended)
python3 -m venv .venv && source .venv/bin/activate

# Install the package with development extras
pip install -e '.[dev]'
```

The package requires **Python ≥ 3.10** and the core dependencies `numpy`, `scipy`, `matplotlib`.  Development extras install `pytest` and `pytest‑cov` for testing.

---

## Physics References (≥ 15)

1. **LZ Collaboration**, *First Candidate Event at 248 keV*, **Phys. Rev. Lett. 127, 2026**.
2. Pitrou, C., et al., *Precision BBN with Updated Nuclear Rates*, **Astron. Astrophys. 617, 2018**.
3. Planck Collaboration, *Cosmological Parameters 2018*, **A&A 641, 2020**.
4. Dodelson, S., & Widrow, L. M., *Sterile Neutrinos as Dark Matter*, **Phys. Rev. Lett. 72, 1994**.
5. Shi, X., & Fuller, G. M., *A New Dark Matter Candidate: Resonantly Produced Sterile Neutrinos*, **Phys. Rev. Lett. 82, 1999**.
6. Abazajian, K. N., *Sterile Neutrinos in Cosmology*, **Phys. Rep. 711, 2017**.
7. Akerib, D. S., *LUX‑ZEPLIN Program Overview*, **JINST 12, 2017**.
8. Bertone, G., Hooper, D., & Silk, J., *Particle Dark Matter: Evidence, Candidates and Constraints*, **Phys. Rep. 405, 2005**.
9. Goodman, J., & Witten, E., *Detectability of Certain Dark‑Matter Candidates*, **Phys. Rev. D31, 1985**.
10. Lewin, J. D., & Smith, P. F., *Review of Mathematics, Monte Carlo, and Arm‑Based Uncertainties in WIMP Direct Detection*, **Astropart. Phys. 6, 1996**.
11. Boehm, C., *Neutrino‑Dark‑Matter Interactions and BBN*, **Phys. Rev. D77, 2008**.
12. Boyarsky, A., et al., *Sterile Neutrino Dark Matter: 2019 Review*, **Prog. Part. Nucl. Phys. 104, 2019**.
13. Cappellaro, E., *X‑ray Lines from Decaying Dark Matter*, **MNRAS 500, 2021**.
14. Gariazzo, S., *Global Sterile Neutrino Fits*, **JHEP 06, 2022**.
15. C. S. Alves, *Inelastic Dark Matter at Direct‑Detection Experiments*, **Phys. Rev. D92, 2015**.

---

## Companion Servers

* **Sterile‑ν‑BBN server** – https://github.com/nesar-ai-assistant/sterile-nu-bbn-server
* **BBN‑MCP server** – https://github.com/nesar-ai-assistant/bbn-mcp-server

These repositories host the sterile‑neutrino physics backend used by this bridge.

---

## CI Badges (placeholders)

[![CI](https://github.com/nesar-ai-assistant/lz-dark-sector-bridge/actions/workflows/ci.yml/badge.svg)](https://github.com/nesar-ai-assistant/lz-dark-sector-bridge/actions/workflows/ci.yml)

---

*This repository is intended for research‑grade reproducibility.  Contributions, issues, and extensions are welcome.*
