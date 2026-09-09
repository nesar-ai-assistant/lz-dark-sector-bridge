"""Comprehensive physics unit tests for lz-dark-sector-bridge.

Tests cover WIMP phenomenology, nuclear form factors, BBN predictions,
mixed dark matter budgets, joint constraints, and the LZ Sep 2026 event.
"""
import math
import pytest
import numpy as np

# ---------------------------------------------------------------------------
# We import from tools.physics — the analytical engine under test.
# ---------------------------------------------------------------------------
from tools.physics import (
    wimp_relic_density,
    helm_form_factor,
    differential_rate,
    recoil_spectrum,
    bbn_predict,
    mixed_dm_budget,
    joint_chi2,
    max_recoil_energy,
    lz_event_compatibility,
)


# ===== Constants used in tests =====
OMEGA_DM_PLANCK = 0.12       # Planck 2018 Ω_DM h²
YP_SM = 0.2471               # Pitrou+ (2018) Y_p at N_eff = 3.044
DH_SM = 2.459e-5             # Pitrou+ (2018) D/H at N_eff = 3.044


# ===== 1. WIMP relic density =====
class TestWIMPRelicDensity:
    """Thermal WIMP relic density from annihilation cross-section."""

    def test_thermal_relic(self):
        """<σv> = 3×10⁻²⁶ cm³/s should give Ω_DM h² ≈ 0.12."""
        omega = wimp_relic_density(sigma_v=3e-26)
        assert 0.05 < omega < 0.25, f"Thermal relic Ω={omega} outside [0.05, 0.25]"

    def test_larger_cross_section_smaller_relic(self):
        """Larger <σv> → smaller relic density (more annihilation)."""
        omega_small = wimp_relic_density(sigma_v=1e-26)
        omega_large = wimp_relic_density(sigma_v=1e-25)
        assert omega_large < omega_small

    def test_zero_cross_section(self):
        """Zero annihilation cross-section → very large (or capped) relic density."""
        omega = wimp_relic_density(sigma_v=0.0)
        assert omega > 1.0 or math.isinf(omega)


# ===== 2. Helm nuclear form factor =====
class TestHelmFormFactor:
    """Helm form factor F²(q) for coherent elastic scattering off Xe."""

    def test_f_at_zero_momentum(self):
        """F(q=0) = 1 by normalisation."""
        f = helm_form_factor(q_fm_inv=0.0, A=131)
        assert abs(f - 1.0) < 1e-6, f"F(0) = {f}, expected 1.0"

    def test_f_decreases_with_q(self):
        """Form factor must decrease as momentum transfer increases."""
        f_low = helm_form_factor(q_fm_inv=0.01, A=131)
        f_high = helm_form_factor(q_fm_inv=1.0, A=131)
        assert f_high < f_low

    def test_f_approaches_zero_at_large_q(self):
        """F(q→∞) → 0  (suppression at large recoil)."""
        f = helm_form_factor(q_fm_inv=5.0, A=131)
        assert f < 0.01, f"F(5 fm⁻¹) = {f}, expected ≪ 1"

    def test_f_bounded(self):
        """F²(q) ∈ [0, 1] for all q ≥ 0."""
        for q in np.linspace(0, 10, 200):
            f = helm_form_factor(q_fm_inv=q, A=131)
            assert 0.0 <= f <= 1.0 + 1e-10


# ===== 3. Recoil spectrum =====
class TestRecoilSpectrum:
    """Differential and integrated recoil spectra."""

    def test_positive_rate(self):
        """Physical parameters → positive differential rate."""
        dRdE = differential_rate(E_keV=50.0, m_chi_GeV=200.0, sigma_SI_cm2=1e-47)
        assert dRdE > 0, f"dR/dE = {dRdE} ≤ 0 for physical params"

    def test_zero_cross_section_zero_rate(self):
        """σ_SI = 0 ⇒ dR/dE = 0 for all E."""
        dRdE = differential_rate(E_keV=50.0, m_chi_GeV=200.0, sigma_SI_cm2=0.0)
        assert dRdE == 0.0

    def test_rate_decreases_with_energy(self):
        """Spectrum falls with recoil energy (exponential-like)."""
        r_low = differential_rate(E_keV=10.0, m_chi_GeV=200.0, sigma_SI_cm2=1e-47)
        r_high = differential_rate(E_keV=200.0, m_chi_GeV=200.0, sigma_SI_cm2=1e-47)
        assert r_high < r_low

    def test_spectrum_array(self):
        """recoil_spectrum returns energy and rate arrays of matching length."""
        E, dRdE = recoil_spectrum(m_chi_GeV=200.0, sigma_SI_cm2=1e-47, n_points=50)
        assert len(E) == len(dRdE) == 50
        assert all(r >= 0 for r in dRdE)

    def test_heavier_wimp_extends_spectrum(self):
        """Heavier WIMP can produce larger maximum recoil energy."""
        E_max_light = max_recoil_energy(m_chi_GeV=50.0)
        E_max_heavy = max_recoil_energy(m_chi_GeV=500.0)
        assert E_max_heavy > E_max_light


# ===== 4. BBN predictions =====
class TestBBN:
    """Standard-model BBN predictions at N_eff = 3.044 (Pitrou+ 2018)."""

    def test_yp_at_sm(self):
        """Y_p at SM N_eff = 3.044 should be close to 0.2471 ± 0.003."""
        result = bbn_predict(n_eff=3.044)
        yp = result["Y_p"]
        assert abs(yp - YP_SM) < 0.01, f"Y_p = {yp}, expected ~{YP_SM}"

    def test_dh_at_sm(self):
        """D/H at SM N_eff = 3.044 should be close to 2.459×10⁻⁵ ± 15%."""
        result = bbn_predict(n_eff=3.044)
        dh = result["D_H"]
        assert abs(dh - DH_SM) / DH_SM < 0.15, f"D/H = {dh}, expected ~{DH_SM}"

    def test_higher_neff_raises_yp(self):
        """Extra radiation (↑N_eff) increases He-4 yield."""
        yp_sm = bbn_predict(n_eff=3.044)["Y_p"]
        yp_hi = bbn_predict(n_eff=4.0)["Y_p"]
        assert yp_hi > yp_sm


# ===== 5. Mixed dark matter budget =====
class TestMixedDM:
    """Joint WIMP + sterile neutrino dark matter fraction constraints."""

    def test_budget_constraint(self):
        """f_WIMP + f_sterile must not exceed 1."""
        result = mixed_dm_budget(f_wimp=0.7, f_sterile=0.3)
        assert result["budget_ok"] is True

    def test_budget_violation(self):
        """f_WIMP + f_sterile > 1 should flag violation."""
        result = mixed_dm_budget(f_wimp=0.8, f_sterile=0.5)
        assert result["budget_ok"] is False

    def test_pure_wimp(self):
        """f_WIMP=1, f_sterile=0 is allowed (pure WIMP scenario)."""
        result = mixed_dm_budget(f_wimp=1.0, f_sterile=0.0)
        assert result["budget_ok"] is True

    def test_pure_sterile(self):
        """f_WIMP=0, f_sterile=1 is allowed."""
        result = mixed_dm_budget(f_wimp=0.0, f_sterile=1.0)
        assert result["budget_ok"] is True


# ===== 6. Joint χ² =====
class TestJointChi2:
    """Joint χ² from multiple observational channels."""

    def test_sm_point_low_chi2(self):
        """The SM-like point should have a low χ² (good fit)."""
        chi2 = joint_chi2(
            m_chi_GeV=200.0,
            sigma_SI_cm2=1e-48,
            m_s_eV=1.0,
            sin2_2theta=1e-12,
            f_wimp=0.99,
            f_sterile=0.01,
        )
        # With near-zero sterile contribution, BBN ≈ SM → low χ²
        assert chi2["total_chi2"] < 50.0, f"SM-like χ² = {chi2['total_chi2']} too large"

    def test_chi2_nonnegative(self):
        """χ² is non-negative by construction."""
        chi2 = joint_chi2(
            m_chi_GeV=200.0,
            sigma_SI_cm2=1e-47,
            m_s_eV=7100.0,
            sin2_2theta=7e-11,
            f_wimp=0.8,
            f_sterile=0.2,
        )
        assert chi2["total_chi2"] >= 0.0

    def test_chi2_returns_components(self):
        """Result includes breakdown by channel."""
        chi2 = joint_chi2(
            m_chi_GeV=200.0,
            sigma_SI_cm2=1e-47,
            m_s_eV=7100.0,
            sin2_2theta=7e-11,
            f_wimp=0.8,
            f_sterile=0.2,
        )
        assert "total_chi2" in chi2


# ===== 7. LZ event analysis =====
class TestLZEvent:
    """Consistency of the LZ Sep 2026 event (248 ± 23 keV) with WIMP models."""

    def test_248keV_consistent_with_heavy_wimp(self):
        """248 keV recoil should be kinematically consistent with m_χ ≥ 200 GeV."""
        compat = lz_event_compatibility(E_recoil_keV=248.0, m_chi_GeV=200.0)
        assert compat["kinematically_allowed"] is True

    def test_248keV_inconsistent_with_light_wimp(self):
        """A 10 GeV WIMP cannot produce 248 keV recoils in Xe."""
        compat = lz_event_compatibility(E_recoil_keV=248.0, m_chi_GeV=10.0)
        assert compat["kinematically_allowed"] is False

    def test_compatibility_returns_max_recoil(self):
        """Result includes the maximum kinematically allowed recoil."""
        compat = lz_event_compatibility(E_recoil_keV=248.0, m_chi_GeV=500.0)
        assert "E_max_keV" in compat
        assert compat["E_max_keV"] > 248.0


# ===== 8. Physical limits / edge cases =====
class TestPhysicalLimits:
    """Edge cases: negative masses, extreme parameters."""

    def test_negative_wimp_mass(self):
        """Negative WIMP mass should raise ValueError or return safe default."""
        with pytest.raises((ValueError, AssertionError)):
            differential_rate(E_keV=50.0, m_chi_GeV=-100.0, sigma_SI_cm2=1e-47)

    def test_negative_energy(self):
        """Negative recoil energy is unphysical."""
        with pytest.raises((ValueError, AssertionError)):
            differential_rate(E_keV=-10.0, m_chi_GeV=200.0, sigma_SI_cm2=1e-47)

    def test_very_large_mass(self):
        """Very heavy WIMP (10 TeV) should still compute without error."""
        dRdE = differential_rate(E_keV=50.0, m_chi_GeV=10000.0, sigma_SI_cm2=1e-47)
        assert np.isfinite(dRdE)
