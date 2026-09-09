"""Tool-level tests for lz-dark-sector-bridge MCP tool wrappers.

These tests exercise the high-level functions in tools/lz_bridge_tools.py
that are registered as MCP tools.  Each function must return an
ArtifactResult dict with at minimum a 'status' key.
"""
import pytest

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


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _is_artifact(result: dict) -> bool:
    """Check that result has the minimal ArtifactResult structure."""
    return isinstance(result, dict) and "status" in result


# ===== describe_tools =====
class TestDescribeTools:
    """The meta-tool that describes all capabilities."""

    def test_returns_ok(self):
        result = describe_lz_bridge_tools()
        assert _is_artifact(result)
        assert result["status"] == "ok"

    def test_lists_tools(self):
        result = describe_lz_bridge_tools()
        assert "tools" in result or "capabilities" in result


# ===== predict_wimp_signal =====
class TestPredictWIMP:
    """WIMP prediction tool wrapper."""

    def test_basic_call(self):
        result = predict_wimp_signal(
            m_chi_GeV=200.0,
            sigma_SI_cm2=1e-47,
            exposure_tonne_yr=5.5,
            E_min_keV=1.5,
            E_max_keV=500.0,
        )
        assert _is_artifact(result)
        assert result["status"] == "ok"

    def test_returns_event_count(self):
        result = predict_wimp_signal(
            m_chi_GeV=200.0,
            sigma_SI_cm2=1e-47,
            exposure_tonne_yr=5.5,
            E_min_keV=1.5,
            E_max_keV=500.0,
        )
        # Must include an expected event count
        assert "expected_events" in result or "n_events" in result

    def test_reasonable_event_count(self):
        """At σ=1e-47 and 5.5 t·yr, expect O(1)–O(10) events."""
        result = predict_wimp_signal(
            m_chi_GeV=200.0,
            sigma_SI_cm2=1e-47,
            exposure_tonne_yr=5.5,
            E_min_keV=1.5,
            E_max_keV=500.0,
        )
        n = result.get("expected_events", result.get("n_events", 0))
        assert 0 < n < 1000, f"Unreasonable event count: {n}"

    def test_zero_exposure_zero_events(self):
        result = predict_wimp_signal(
            m_chi_GeV=200.0,
            sigma_SI_cm2=1e-47,
            exposure_tonne_yr=0.0,
            E_min_keV=1.5,
            E_max_keV=500.0,
        )
        n = result.get("expected_events", result.get("n_events", 0))
        assert n == 0.0 or n == pytest.approx(0.0, abs=1e-10)


# ===== predict_mixed_dm =====
class TestPredictMixedDM:
    """Mixed dark matter budget tool."""

    def test_valid_budget(self):
        result = predict_mixed_dm(
            f_wimp=0.7,
            f_sterile=0.3,
            m_chi_GeV=200.0,
            sigma_SI_cm2=1e-47,
            m_s_eV=7100.0,
            sin2_2theta=7e-11,
        )
        assert _is_artifact(result)
        assert result["status"] == "ok"

    def test_budget_violation_flagged(self):
        result = predict_mixed_dm(
            f_wimp=0.8,
            f_sterile=0.5,
            m_chi_GeV=200.0,
            sigma_SI_cm2=1e-47,
            m_s_eV=7100.0,
            sin2_2theta=7e-11,
        )
        assert _is_artifact(result)
        # Should flag the violation somewhere
        assert result.get("budget_ok") is False or "warning" in str(result).lower()


# ===== predict_recoil_spectrum =====
class TestPredictRecoilSpectrum:
    """Recoil spectrum tool."""

    def test_returns_spectrum_data(self):
        result = predict_recoil_spectrum(
            m_chi_GeV=200.0,
            sigma_SI_cm2=1e-47,
            n_points=20,
            E_min_keV=1.5,
            E_max_keV=500.0,
        )
        assert _is_artifact(result)
        assert result["status"] == "ok"


# ===== compute_joint_constraints =====
class TestJointConstraints:
    """Joint constraint χ² tool."""

    def test_returns_chi2(self):
        result = compute_joint_constraints(
            m_chi_GeV=200.0,
            sigma_SI_cm2=1e-47,
            m_s_eV=7100.0,
            sin2_2theta=7e-11,
            f_wimp=0.8,
            f_sterile=0.2,
        )
        assert _is_artifact(result)
        assert result["status"] == "ok"
        assert "total_chi2" in result or "chi2" in result


# ===== evaluate_lz_event =====
class TestEvaluateLZEvent:
    """LZ Sep 2026 event analysis tool."""

    def test_returns_analysis(self):
        result = evaluate_lz_event(
            E_recoil_keV=248.0,
            E_sigma_keV=23.0,
            m_chi_min_GeV=50.0,
            m_chi_max_GeV=5000.0,
            n_mass_points=10,
        )
        assert _is_artifact(result)
        assert result["status"] == "ok"


# ===== plot tools (may create files) =====
class TestPlotTools:
    """Plot tools should return paths to generated files."""

    def test_plot_exclusion(self, tmp_path):
        result = plot_wimp_exclusion(
            m_min_GeV=10.0,
            m_max_GeV=1e4,
            n_points=10,
            output_dir=str(tmp_path),
        )
        assert _is_artifact(result)

    def test_plot_mixed_dm(self, tmp_path):
        result = plot_mixed_dm_plane(
            n_grid=5,
            output_dir=str(tmp_path),
        )
        assert _is_artifact(result)
