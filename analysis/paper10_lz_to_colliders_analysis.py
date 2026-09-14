#!/usr/bin/env python3
"""
Paper 10 Analysis: From LUX‑Z E to Colliders – Higgsino Dark Matter
S. Di Mauro & J. Shaikh (2026), arXiv:2609.08712
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add the repo root to sys.path for importing physics utilities
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tools.physics import (
    wimp_relic_density,
    helm_form_factor,
    OBS_DATA,
)

# Output directory for generated figures
OUTDIR = os.path.join(os.path.dirname(__file__), 'output')
os.makedirs(OUTDIR, exist_ok=True)

# ------------------------------------------------------------
# 1. Higgsino relic density vs mass (thermal freeze‑out)
# ------------------------------------------------------------

masses = np.linspace(800, 1500, 200)  # GeV
# Approximate annihilation cross‑section for pure Higgsino (dominant to WW)
# sigma_v ≈ 2e-26 cm^3/s * (1 TeV / m)^2  (scales roughly as 1/m^2)
sigma_vs = 2e-26 * (1000.0 / masses) ** 2
relics = np.array([wimp_relic_density(s) for s in sigma_vs])

plt.figure(figsize=(8, 5))
plt.semilogy(masses, relics, 'C0', lw=2)
plt.axhline(OBS_DATA['Omega_DM'].value, color='gold', ls='--', lw=1,
            label='Planck Ω_DM h^2')
plt.xlabel('Higgsino mass $m_{\chi}$ [GeV]')
plt.ylabel('Thermal relic density $\Omega_\chi h^2$')
plt.title('Thermal relic density of a pure Higgsino')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'higgsino_relic_vs_mass.png'))
plt.close()

# ------------------------------------------------------------
# 2. Collider reach: chargino‑neutralino mass splitting vs LHC limit
# ------------------------------------------------------------
# Simplified: chargino–neutralino splitting Δm ≈ 350 MeV for pure Higgsino.
# LHC excludes m_χ < 1150 GeV for Δm ≈ 350 MeV (extrapolated from ATLAS/ CMS HL‑LHC studies).

m_vals = np.linspace(800, 1500, 200)
# Assume a step function exclusion: excluded if m < 1150 GeV
excl = m_vals < 1150

plt.figure(figsize=(8, 5))
plt.plot(m_vals, excl.astype(int), drawstyle='steps-post', color='C3')
plt.xlabel('Higgsino mass $m_{\chi}$ [GeV]')
plt.ylabel('Collider exclusion (1 = excluded)')
plt.title('Approximate LHC exclusion for pure Higgsino (Δm≈350 MeV)')
plt.ylim(-0.1, 1.1)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'higgsino_collider_exclusion.png'))
plt.close()

print('Analysis complete. Figures saved to', OUTDIR)
