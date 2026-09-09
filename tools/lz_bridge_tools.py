"""LZ Bridge tool wrappers for the MCP server.

Each function mirrors the low‑level physics functions from ``tools.physics``
but packages the result into a minimal *ArtifactResult* dict with a ``status``
field (required by the test suite).  The implementation is deliberately
lightweight – the heavy physics lives in ``tools.physics`` – but the
dictionary structure matches what the client side expects when calling the
MCP server via ``describe_tools``.
"""

from __future__ import annotations

from typing import Dict, Any, Tuple

import json
import numpy as np

from .physics import (
    wimp_relic_density,
    helm_form_factor,
    differential_rate,
    recoil_spectrum,
    bbn_predict,
    mixed_dm_budget,
    joint_chi2,
    lz_event_compatibility,
    expected_lz_events,
    wimp_delta_neff,
    sterile_delta_neff,
    sterile_relic_density,
    radiative_decay_rate,
    xray_flux,
    free_streaming_length_kpc,
    chi2_bbn,
)

# ---------------------------------------------------------------------------
# Helper to format a successful ArtifactResult
# ---------------------------------------------------------------------------
def _artifact(status: str = "ok", files: list = None, message: str = "", metadata: dict = None) -> dict:
    return {
        "status": status,
        "files": files or [],
        "message": message,
        "metadata": metadata or {},
    }

# ---------------------------------------------------------------------------
# 1. describe_lz_bridge_tools – meta information
# ---------------------------------------------------------------------------
def describe_lz_bridge_tools() -> dict:
    """Return a description of the bridge toolset.

    Mirrors the ArtifactResult pattern used by the sterile‑nu‑bbn server.
    """
    tools = [
        {"name": "describe_lz_bridge_tools", "signature": "() -> dict"},
        {"name": "predict_wimp_signal", "signature": "(m_chi_GeV, sigma_SI_cm2, exposure_tonne_yr=2.84, E_min_keV=5.0, E_max_keV=270.0) -> dict"},
        {"name": "predict_mixed_dm", "signature": "(f_wimp, f_sterile, m_chi_GeV, sigma_SI_cm2, m_s_eV, sin2_2theta) -> dict"},
        {"name": "predict_recoil_spectrum", "signature": "(m_chi_GeV, sigma_SI_cm2, n_points=100, E_min_keV=1.5, E_max_keV=500.0) -> dict"},
        {"name": "compute_joint_constraints", "signature": "(m_chi_GeV, sigma_SI_cm2, m_s_eV, sin2_2theta, f_wimp, f_sterile) -> dict"},
        {"name": "evaluate_lz_event", "signature": "(E_recoil_keV=248.0, E_sigma_keV=23.0, m_chi_min_GeV=50.0, m_chi_max_GeV=5000.0, n_mass_points=50) -> dict"},
        {"name": "plot_wimp_exclusion", "signature": "(m_min_GeV=10, m_max_GeV=1e4, n_points=80, output_dir='output') -> dict"},
        {"name": "plot_mixed_dm_plane", "signature": "(n_grid=40, output_dir='output') -> dict"},
    ]
    message = "LZ Dark Sector Bridge – toolset description."
    metadata = {"tools": tools, "description": "Bridge between WIMP and sterile‑neutrino physics for LZ."}
    # Also expose top‑level key for test compatibility
    result = _artifact(message=message, metadata=metadata)
    result["tools"] = tools
    return result

# ---------------------------------------------------------------------------
# 2. predict_wimp_signal – event count and spectrum
# ---------------------------------------------------------------------------
def predict_wimp_signal(
    m_chi_GeV: float,
    sigma_SI_cm2: float,
    exposure_tonne_yr: float = 2.84,
    E_min_keV: float = 5.0,
    E_max_keV: float = 270.0,
) -> dict:
    """Calculate expected LZ events and return recoil spectrum.

    Returns an ArtifactResult with ``expected_events`` top‑level for test
    compatibility and stores full data in ``metadata``.
    """
    # Expected events using physics helper
    expected = expected_lz_events(
        m_chi_GeV, sigma_SI_cm2, E_min_keV, E_max_keV
    )
    # Scale by exposure ratio if different from default 2.84 tonne-yr
    expected *= exposure_tonne_yr / 2.84
    # Full recoil spectrum (for possible downstream use)
    Es, rates = recoil_spectrum(m_chi_GeV, sigma_SI_cm2, n_points=200, E_min_keV=E_min_keV, E_max_keV=E_max_keV)
    metadata = {
        "expected_events": expected,
        "energies_keV": Es.tolist(),
        "rates_per_keV": rates.tolist(),
        "wimp_delta_neff": wimp_delta_neff(m_chi_GeV),
    }
    result = _artifact(message="WIMP signal prediction", metadata=metadata)
    # expose expected_events at top level for the test suite
    result["expected_events"] = expected
    return result

# ---------------------------------------------------------------------------
# 3. predict_mixed_dm – budget evaluation and relic densities
# ---------------------------------------------------------------------------
def predict_mixed_dm(
    f_wimp: float,
    f_sterile: float,
    m_chi_GeV: float,
    sigma_SI_cm2: float,
    m_s_eV: float,
    sin2_2theta: float,
) -> dict:
    """Assess mixed DM scenario and return budget flag and relic densities.

    The result follows the ArtifactResult format with top‑level ``budget_ok``
    for easy test access.
    """
    budget = mixed_dm_budget(f_wimp, f_sterile)
    wimp_relic = wimp_relic_density(sigma_SI_cm2)
    sterile_relic = sterile_relic_density(m_s_eV, sin2_2theta)
    metadata = {
        "budget_ok": budget["budget_ok"],
        "wimp_relic_density": wimp_relic,
        "sterile_relic_density": sterile_relic,
        "f_wimp": f_wimp,
        "f_sterile": f_sterile,
    }
    result = _artifact(message="Mixed DM budget", metadata=metadata)
    result["budget_ok"] = budget["budget_ok"]
    return result

# ---------------------------------------------------------------------------
# 4. predict_recoil_spectrum – thin wrapper
# ---------------------------------------------------------------------------
def predict_recoil_spectrum(
    m_chi_GeV: float,
    sigma_SI_cm2: float,
    n_points: int = 100,
    E_min_keV: float = 1.5,
    E_max_keV: float = 500.0,
) -> dict:
    Es, rates = recoil_spectrum(m_chi_GeV, sigma_SI_cm2, n_points, E_min_keV, E_max_keV)
    metadata = {
        "energies_keV": Es.tolist(),
        "rates": rates.tolist(),
    }
    return _artifact(message="Recoil spectrum", metadata=metadata)

# ---------------------------------------------------------------------------
# 5. compute_joint_constraints – summary chi2
# ---------------------------------------------------------------------------
def compute_joint_constraints(
    m_chi_GeV: float,
    sigma_SI_cm2: float,
    m_s_eV: float,
    sin2_2theta: float,
    f_wimp: float,
    f_sterile: float,
) -> dict:
    chi2_dict = joint_chi2(m_chi_GeV, sigma_SI_cm2, m_s_eV, sin2_2theta, f_wimp, f_sterile)
    # expose total chi2 at top level for test expectations
    result = _artifact(message="Joint chi2", metadata=chi2_dict)
    if "total_chi2" in chi2_dict:
        result["total_chi2"] = chi2_dict["total_chi2"]
    return result

# ---------------------------------------------------------------------------
# 6. evaluate_lz_event – scan WIMP masses for compatibility
# ---------------------------------------------------------------------------
def evaluate_lz_event(
    E_recoil_keV: float = 248.0,
    E_sigma_keV: float = 23.0,
    m_chi_min_GeV: float = 50.0,
    m_chi_max_GeV: float = 5000.0,
    n_mass_points: int = 50,
) -> dict:
    masses = np.linspace(m_chi_min_GeV, m_chi_max_GeV, n_mass_points)
    allowed = []
    for m in masses:
        compat = lz_event_compatibility(E_recoil_keV, m)
        allowed.append(compat["kinematically_allowed"])
    metadata = {
        "masses_GeV": masses.tolist(),
        "allowed": allowed,
        "E_recoil_keV": E_recoil_keV,
        "E_sigma_keV": E_sigma_keV,
    }
    return _artifact(message="LZ event compatibility scan", metadata=metadata)

# ---------------------------------------------------------------------------
# 7. Plot helpers – generate real exclusion plot (300 DPI)
# ---------------------------------------------------------------------------
def _ensure_dir(path: str) -> str:
    import os
    os.makedirs(path, exist_ok=True)
    return path

def plot_wimp_exclusion(
    m_min_GeV: float = 10.0,
    m_max_GeV: float = 1e4,
    n_points: int = 80,
    output_dir: str = "output",
) -> dict:
    """Create a publication‑quality exclusion curve for WIMPs.

    The exclusion is defined heuristically: points where the expected number of
    events exceeds three are considered excluded.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = _ensure_dir(output_dir)
    masses = np.logspace(np.log10(m_min_GeV), np.log10(m_max_GeV), n_points)
    # For each mass we find the cross‑section that yields ~3 events using a simple scaling.
    sigma_ref = 1e-45
    sigmas = []
    for m in masses:
        ev_ref = expected_lz_events(m, sigma_ref)
        if ev_ref == 0:
            sig = float('inf')
        else:
            sig = sigma_ref * (3.0 / ev_ref)
        sigmas.append(sig)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.loglog(masses, sigmas, label="3‑event exclusion (approx)")
    ax.set_xlabel(r"WIMP mass $m_\chi$ [GeV]")
    ax.set_ylabel(r"Spin‑independent cross‑section $\sigma_{SI}$ [cm$^2$]")
    ax.set_title("LZ projected 3‑event exclusion curve")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend()
    fig_path = f"{out_dir}/wimp_exclusion.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return _artifact(message="WIMP exclusion plot", files=[fig_path], metadata={"plot_path": fig_path})

# ---------------------------------------------------------------------------
# 8. Mixed‑DM plane plot – chi2 map
# ---------------------------------------------------------------------------
def plot_mixed_dm_plane(
    n_grid: int = 40,
    output_dir: str = "output",
) -> dict:
    """Generate a contour plot of the joint chi² in the (f_wimp, f_sterile) plane.

    The grid spans f_wimp and f_sterile from 0 to 1.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = _ensure_dir(output_dir)
    f_vals = np.linspace(0, 1, n_grid)
    Fw, Fs = np.meshgrid(f_vals, f_vals)
    chi2_grid = np.empty_like(Fw)
    # Representative parameters for the other variables
    m_chi, sigma, m_s, theta = 200.0, 1e-47, 7.1e3, 7e-11
    for i in range(n_grid):
        for j in range(n_grid):
            chi = joint_chi2(m_chi, sigma, m_s, theta, float(Fw[i, j]), float(Fs[i, j]))
            chi2_grid[i, j] = chi.get("total_chi2", np.nan)
    fig, ax = plt.subplots(figsize=(7, 6))
    cs = ax.contourf(Fw, Fs, chi2_grid, levels=[0, 5, 10, 20, 50], cmap="viridis")
    fig.colorbar(cs, label=r"Total $\chi^2$")
    ax.set_xlabel("f_wimp")
    ax.set_ylabel("f_sterile")
    ax.set_title(r"Joint $\chi^2$ in mixed DM fraction space")
    fig_path = f"{out_dir}/mixed_dm_plane.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return _artifact(message="Mixed DM plane plot", files=[fig_path], metadata={"plot_path": fig_path})

# End of file