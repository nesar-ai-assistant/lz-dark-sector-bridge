"""MCP server entry point for lz-dark-sector-bridge (MCP SDK v2)."""
from __future__ import annotations

from mcp.server.mcpserver import MCPServer

from tools.lz_bridge_tools import (
    describe_lz_bridge_tools,
    predict_wimp_signal,
    predict_mixed_dm,
    predict_recoil_spectrum,
    compute_joint_constraints,
    evaluate_lz_event,
    plot_wimp_exclusion,
    plot_mixed_dm_plane,
)

mcp = MCPServer(
    name="lz-dark-sector-bridge",
    description=(
        "Bridges LZ direct dark-matter detection with sterile neutrino BBN "
        "cosmology.  Predicts WIMP recoil spectra, mixed dark-matter budgets, "
        "joint χ² constraints, and compares against the LZ Sep 2026 result."
    ),
)


@mcp.tool()
def describe_tools() -> str:
    """Describe this server's capabilities, the physics, and the observational data it uses. Call this first."""
    import json
    result = describe_lz_bridge_tools()
    return json.dumps(result, indent=2, default=str)


@mcp.tool()
def predict_wimp(
    m_chi_GeV: float = 200.0,
    sigma_SI_cm2: float = 1e-47,
    exposure_tonne_yr: float = 5.5,
    E_min_keV: float = 1.5,
    E_max_keV: float = 500.0,
) -> str:
    """Predict WIMP direct-detection signal for given mass and cross-section.

    Parameters:
        m_chi_GeV: WIMP mass in GeV (default: 200)
        sigma_SI_cm2: spin-independent cross-section in cm² (default: 1e-47)
        exposure_tonne_yr: detector exposure in tonne·yr (default: 5.5, LZ full)
        E_min_keV: recoil energy threshold in keV (default: 1.5)
        E_max_keV: maximum recoil energy in keV (default: 500)
    """
    import json
    result = predict_wimp_signal(m_chi_GeV, sigma_SI_cm2, exposure_tonne_yr, E_min_keV, E_max_keV)
    return json.dumps(result, indent=2, default=str)


@mcp.tool()
def predict_mixed_dm_budget(
    f_wimp: float = 0.8,
    f_sterile: float = 0.2,
    m_chi_GeV: float = 200.0,
    sigma_SI_cm2: float = 1e-47,
    m_s_eV: float = 7100.0,
    sin2_2theta: float = 7e-11,
) -> str:
    """Evaluate a mixed WIMP + sterile neutrino dark-matter budget.

    Parameters:
        f_wimp: WIMP fraction of Ω_DM (default: 0.8)
        f_sterile: sterile neutrino fraction of Ω_DM (default: 0.2)
        m_chi_GeV: WIMP mass in GeV
        sigma_SI_cm2: SI cross-section in cm²
        m_s_eV: sterile neutrino mass in eV
        sin2_2theta: sterile mixing parameter sin²(2θ)
    """
    import json
    result = predict_mixed_dm(f_wimp, f_sterile, m_chi_GeV, sigma_SI_cm2, m_s_eV, sin2_2theta)
    return json.dumps(result, indent=2, default=str)


@mcp.tool()
def recoil_spectrum(
    m_chi_GeV: float = 200.0,
    sigma_SI_cm2: float = 1e-47,
    n_points: int = 100,
    E_min_keV: float = 1.5,
    E_max_keV: float = 500.0,
) -> str:
    """Compute the differential recoil energy spectrum dR/dE for xenon.

    Parameters:
        m_chi_GeV: WIMP mass in GeV
        sigma_SI_cm2: SI cross-section in cm²
        n_points: number of energy bins
        E_min_keV: minimum recoil energy in keV
        E_max_keV: maximum recoil energy in keV
    """
    import json
    result = predict_recoil_spectrum(m_chi_GeV, sigma_SI_cm2, n_points, E_min_keV, E_max_keV)
    return json.dumps(result, indent=2, default=str)


@mcp.tool()
def joint_constraints(
    m_chi_GeV: float = 200.0,
    sigma_SI_cm2: float = 1e-47,
    m_s_eV: float = 7100.0,
    sin2_2theta: float = 7e-11,
    f_wimp: float = 0.8,
    f_sterile: float = 0.2,
) -> str:
    """Compute joint χ² from LZ + BBN + Planck + X-ray constraints.

    Parameters:
        m_chi_GeV: WIMP mass in GeV
        sigma_SI_cm2: SI cross-section in cm²
        m_s_eV: sterile neutrino mass in eV
        sin2_2theta: sterile mixing parameter
        f_wimp: WIMP fraction of Ω_DM
        f_sterile: sterile neutrino fraction
    """
    import json
    result = compute_joint_constraints(
        m_chi_GeV, sigma_SI_cm2, m_s_eV, sin2_2theta, f_wimp, f_sterile
    )
    return json.dumps(result, indent=2, default=str)


@mcp.tool()
def lz_event_analysis(
    E_recoil_keV: float = 248.0,
    E_sigma_keV: float = 23.0,
    m_chi_min_GeV: float = 50.0,
    m_chi_max_GeV: float = 5000.0,
    n_mass_points: int = 50,
) -> str:
    """Analyse the LZ Sep 2026 event: which WIMP masses are consistent?

    Parameters:
        E_recoil_keV: observed recoil energy (default: 248 keV)
        E_sigma_keV: energy uncertainty (default: 23 keV)
        m_chi_min_GeV: minimum WIMP mass to scan
        m_chi_max_GeV: maximum WIMP mass to scan
        n_mass_points: number of mass points
    """
    import json
    result = evaluate_lz_event(E_recoil_keV, E_sigma_keV, m_chi_min_GeV, m_chi_max_GeV, n_mass_points)
    return json.dumps(result, indent=2, default=str)


@mcp.tool()
def plot_wimp_exclusion_curve(
    m_min_GeV: float = 10.0,
    m_max_GeV: float = 1e4,
    n_points: int = 80,
) -> str:
    """Plot the WIMP mass vs cross-section exclusion curve. Returns path to PNG.

    Parameters:
        m_min_GeV: minimum WIMP mass
        m_max_GeV: maximum WIMP mass
        n_points: number of mass points
    """
    import json
    result = plot_wimp_exclusion(m_min_GeV, m_max_GeV, n_points, output_dir="output")
    return json.dumps(result, indent=2, default=str)


@mcp.tool()
def plot_mixed_dm(
    n_grid: int = 40,
) -> str:
    """Plot the (f_WIMP, f_sterile) plane with constraint regions. Returns path to PNG.

    Parameters:
        n_grid: grid resolution per axis
    """
    import json
    result = plot_mixed_dm_plane(n_grid, output_dir="output")
    return json.dumps(result, indent=2, default=str)


if __name__ == "__main__":
    mcp.run(transport="stdio")
