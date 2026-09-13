#!/usr/bin/env python3
"""
Paper 9 Analysis: Solar Capture Tests of Inelastic Dark Matter after LZ
Di Mauro & Shaikh (2026), arXiv:2609.06760
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from tools.physics import (
    bbn_predict, sterile_delta_neff, wimp_delta_neff,
    joint_chi2, mixed_dm_budget
)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 12, 'axes.labelsize': 14, 'axes.titlesize': 15,
    'legend.fontsize': 10, 'xtick.labelsize': 11, 'ytick.labelsize': 11,
    'figure.dpi': 100, 'savefig.dpi': 300, 'figure.facecolor': 'white',
    'axes.grid': True, 'grid.alpha': 0.3,
})

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
os.makedirs(OUTDIR, exist_ok=True)

t_sun = 4.57e9 * 3.156e7
R_sun = 6.957e10
sigma_v_hig = 1.3e-26
Gamma_lim_WW = 1.5e19
Gamma_lim_bb = 7.5e19
Gamma_lim_WW_cons = 2.44e19

def solar_capture_model(delta_keV, C0, delta_LZ, delta_max, delta_scale):
    delta = np.atleast_1d(np.float64(delta_keV))
    raw = np.where(delta < delta_max,
        (1.0 - delta / delta_max)**2 * np.exp(-delta / delta_scale), 0.0)
    raw_at_LZ = (1.0 - delta_LZ / delta_max)**2 * np.exp(-delta_LZ / delta_scale)
    norm = C0 / raw_at_LZ if raw_at_LZ > 0 else 0.0
    return norm * raw

def annihilation_rate(C_sun, R_chi_cm, sigma_v=sigma_v_hig):
    V_eff = (4.0 / 3.0) * np.pi * R_chi_cm**3
    C_A = sigma_v / V_eff
    C_sun = np.atleast_1d(np.float64(C_sun))
    arg = C_sun * C_A
    tau_eq = np.where(arg > 0, 1.0 / np.sqrt(arg), 1e30)
    ratio = t_sun / tau_eq
    tanh_val = np.where(ratio > 50, 1.0, np.tanh(ratio))
    return 0.5 * C_sun * tanh_val**2

# ===================== PLOT 1 =====================
print("Plot 1: Solar capture rate vs splitting...")
delta_arr = np.linspace(200, 600, 500)
C_higgsino = solar_capture_model(delta_arr, 1.02e23, 377, 580, 120)
C_pseudo = solar_capture_model(delta_arr, 1.48e20, 297, 450, 100)
C_neutron = solar_capture_model(delta_arr, 3.98e17, 300, 420, 90)

fig, ax = plt.subplots(figsize=(10, 7))
ax.semilogy(delta_arr, C_higgsino, 'C0', lw=2.5,
    label=r'Thermal Higgsino ($m_\chi=1.08$ TeV)')
ax.semilogy(delta_arr, C_pseudo, 'C1', lw=2.5,
    label=r'Pseudo-Dirac vector ($m_\chi=1$ TeV)')
ax.semilogy(delta_arr, C_neutron, 'C2', lw=2.5,
    label=r'Neutron-philic $\mathcal{O}_4$ ($m_\chi=1$ TeV)')
ax.axhline(2*Gamma_lim_WW, color='cyan', lw=2, ls='-', alpha=0.9,
    label=r'$2\Gamma_{\rm lim}^{WW}$ (IceCube)')
ax.axhline(2*Gamma_lim_bb, color='cyan', lw=2, ls='--', alpha=0.7,
    label=r'$2\Gamma_{\rm lim}^{b\bar{b}}$ (IceCube)')
for dval, col in [(377, 'C0'), (297, 'C1'), (300, 'C2')]:
    ax.axvline(dval, color=col, ls=':', lw=1.5, alpha=0.6)
ax.axvspan(250, 400, alpha=0.08, color='gold', label='LZ preferred region')
ax.plot(377, 1.02e23, 'o', color='C0', ms=10, mec='k', zorder=5)
ax.plot(297, 1.48e20, 's', color='C1', ms=10, mec='k', zorder=5)
ax.plot(300, 3.98e17, 'D', color='C2', ms=10, mec='k', zorder=5)
ax.set_xlabel(r'Endothermic splitting $\delta$ [keV]')
ax.set_ylabel(r'Solar capture rate $C_\odot$ [s$^{-1}$]')
ax.set_title('Solar Capture Rate vs. Endothermic Splitting\n(Di Mauro & Shaikh 2026, arXiv:2609.06760)')
ax.set_xlim(200, 600); ax.set_ylim(1e10, 1e26)
ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, 'solar_capture_vs_splitting.png'))
plt.close(fig)
print("  -> solar_capture_vs_splitting.png")

# ===================== PLOT 2 =====================
print("Plot 2: Higgsino exclusion...")
delta_arr2 = np.linspace(300, 590, 500)
C_hig = solar_capture_model(delta_arr2, 1.02e23, 377, 580, 120)
R_stall = 0.200 * R_sun
R_SD = 0.0228 * R_sun
R_therm = 0.0037 * R_sun
Gamma_stall = annihilation_rate(C_hig, R_stall)
Gamma_SD = annihilation_rate(C_hig, R_SD)
Gamma_therm = annihilation_rate(C_hig, R_therm)

fig, ax = plt.subplots(figsize=(10, 7))
ax.semilogy(delta_arr2, Gamma_stall, 'k', lw=2.5,
    label=r'Tree-level stall ($R_\chi \approx 0.20\,R_\odot$)')
ax.semilogy(delta_arr2, Gamma_SD, 'C0', lw=2.5, ls='--',
    label=r'SD loop cooling ($R_\chi \approx 0.023\,R_\odot$)')
ax.semilogy(delta_arr2, Gamma_therm, 'C3', lw=2.5, ls='-.',
    label=r'SI+SD thermalized ($R_\chi \approx 0.0037\,R_\odot$)')
ax.axhline(Gamma_lim_WW, color='cyan', lw=2.5, alpha=0.9,
    label=r'IceCube $WW$ limit')
ax.axhline(Gamma_lim_WW_cons, color='cyan', lw=1.5, ls='--', alpha=0.6,
    label=r'$WW$-only mapping')
ax.fill_between(delta_arr2, Gamma_lim_WW, 1e23, alpha=0.10, color='red')
ax.text(310, 3e19, 'EXCLUDED by IceCube', fontsize=13, fontweight='bold',
    color='red', alpha=0.7)
ax.axvline(377, color='C0', ls=':', lw=2, alpha=0.7)
ax.text(380, 2e11, r'$\delta_{\rm LZ} \approx 377$ keV', fontsize=10,
    color='C0', rotation=90, va='bottom')
for dl, lab, col in [(512.5, r'$\delta_{\rm lim}^{\rm stall}$', 'k'),
    (558.6, r'$\delta_{\rm lim}^{\rm SD}$', 'C0'),
    (564.8, r'$\delta_{\rm lim}^{\rm therm}$', 'C3')]:
    ax.axvline(dl, color=col, ls=':', lw=1.2, alpha=0.5)
    ax.text(dl+2, 5e10, lab, fontsize=9, color=col, rotation=90, va='bottom')
ax.annotate('', xy=(377, 5e17), xytext=(512, 5e17),
    arrowprops=dict(arrowstyle='<->', color='darkred', lw=2))
ax.text(445, 8e17, r'$\delta_{\rm LZ} \ll \delta_{\rm lim}$'+'\nHiggsino excluded',
    fontsize=10, ha='center', color='darkred', fontweight='bold')
ax.set_xlabel(r'Endothermic splitting $\delta$ [keV]')
ax.set_ylabel(r'Annihilation rate $\Gamma_A(t_\odot)$ [s$^{-1}$]')
ax.set_title('Higgsino Solar Neutrino Constraint\n(Three post-capture cooling scenarios)')
ax.set_xlim(300, 590); ax.set_ylim(1e10, 5e22)
ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, 'higgsino_exclusion_delta.png'))
plt.close(fig)
print("  -> higgsino_exclusion_delta.png")

# ===================== PLOT 3 =====================
print("Plot 3: Mixed DM constraints...")
n_grid = 60
f_wimp_arr = np.linspace(0.01, 1.0, n_grid)
f_sterile_arr = np.linspace(0.01, 1.0, n_grid)
FW, FS = np.meshgrid(f_wimp_arr, f_sterile_arr)
m_chi, sigma_SI, m_s, sin2_2theta = 1080.0, 4e-50, 7000.0, 7e-11
chi2_grid = np.zeros_like(FW)
for i in range(n_grid):
    for j in range(n_grid):
        try:
            result = joint_chi2(m_chi, sigma_SI, m_s, sin2_2theta, FW[i,j], FS[i,j])
            chi2_grid[i,j] = result['total_chi2']
        except:
            chi2_grid[i,j] = np.nan
chi2_grid = np.clip(chi2_grid, 0.01, 200)

fig, ax = plt.subplots(figsize=(9, 7))
pcm = ax.pcolormesh(FW, FS, chi2_grid, cmap='RdYlGn_r', shading='auto', vmin=0, vmax=50)
fig.colorbar(pcm, ax=ax, label=r'Total $\chi^2$', shrink=0.85)
fw_line = np.linspace(0, 1, 100)
ax.plot(fw_line, 1.0-fw_line, 'k--', lw=2, label=r'$f_{\rm WIMP}+f_{\rm sterile}=1$')
ax.axvline(1.0, color='red', lw=3, label='Pure Higgsino (IceCube excluded)')
ax.fill_betweenx([0,1], 0.85, 1.0, alpha=0.15, color='red')
ax.text(0.88, 0.5, 'Solar\ncapture\nexcluded', fontsize=10, color='red',
    ha='center', va='center', fontweight='bold', rotation=90)
ax.axvspan(0.12, 0.15, alpha=0.2, color='blue',
    label=r'PQ subcomponent ($f_\chi \approx 0.12$--$0.15$)')
ax.axhspan(0.8, 1.0, alpha=0.1, color='purple')
ax.text(0.5, 0.9, r'Ly-$\alpha$ tension for 7 keV sterile', fontsize=9,
    ha='center', color='purple', alpha=0.7)
ax.annotate(r'Preferred: $f_{\rm WIMP}\sim 0.1$--$0.3$', xy=(0.2, 0.5),
    fontsize=11, fontweight='bold',
    bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='gray', alpha=0.8))
ax.set_xlabel(r'WIMP fraction $f_{\rm WIMP}$')
ax.set_ylabel(r'Sterile neutrino fraction $f_{\rm sterile}$')
ax.set_title('Mixed DM Constraints After Higgsino Solar Exclusion\n'
    r'($m_\chi=1.08$ TeV, $m_s=7$ keV, $\sin^2 2\theta=7\times10^{-11}$)')
ax.set_xlim(0, 1.05); ax.set_ylim(0, 1.05)
ax.legend(loc='lower left', fontsize=8.5, framealpha=0.9)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, 'mixed_dm_after_higgsino_exclusion.png'))
plt.close(fig)
print("  -> mixed_dm_after_higgsino_exclusion.png")

# ===================== PLOT 4 =====================
print("Plot 4: BBN predictions...")
Neff_arr = np.linspace(2.5, 4.5, 200)
Yp_arr = np.array([bbn_predict(n)['Y_p'] for n in Neff_arr])
DH_arr = np.array([bbn_predict(n)['D_H'] for n in Neff_arr])
dN_sterile = sterile_delta_neff(7000, 7e-11)
dN_wimp = wimp_delta_neff(1080)
Neff_SM = 3.044
Neff_s = Neff_SM + dN_sterile
Neff_ws = Neff_SM + dN_wimp + dN_sterile
Yp_obs, Yp_err = 0.2449, 0.0040
DH_obs, DH_err = 2.547e-5, 0.025e-5
Neff_planck, Neff_planck_err = 2.99, 0.17

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
ax1.plot(Neff_arr, Yp_arr, 'k-', lw=2.5, label='BBN prediction (Pitrou+ 2018)')
ax1.axhspan(Yp_obs-Yp_err, Yp_obs+Yp_err, alpha=0.25, color='gold',
    label=r'$Y_p = 0.2449 \pm 0.0040$')
ax1.axhline(Yp_obs, color='goldenrod', ls='--', lw=1)
ax1.axvspan(Neff_planck-2*Neff_planck_err, Neff_planck+2*Neff_planck_err,
    alpha=0.12, color='skyblue', label=r'Planck $N_{\rm eff}$ ($2\sigma$)')
ax1.axvline(Neff_SM, color='green', lw=2, alpha=0.7, label=r'SM: $N_{\rm eff}=3.044$')
ax1.axvline(Neff_s, color='blue', lw=2, ls='--', alpha=0.7,
    label=r'SM + 7 keV $\nu_s$: $N_{\rm eff}=%.3f$' % Neff_s)
ax1.axvline(Neff_ws, color='red', lw=2, ls=':', alpha=0.7,
    label=r'SM + WIMP + $\nu_s$: $N_{\rm eff}=%.3f$' % Neff_ws)
ax1.set_ylabel(r'$Y_p$ (primordial $^4$He)')
ax1.set_title('BBN Predictions for LZ-Motivated Dark Matter Scenarios', fontsize=14)
ax1.legend(loc='upper left', fontsize=8.5, framealpha=0.9)
ax1.set_ylim(0.235, 0.260)

ax2.plot(Neff_arr, DH_arr*1e5, 'k-', lw=2.5, label='BBN prediction')
ax2.axhspan((DH_obs-DH_err)*1e5, (DH_obs+DH_err)*1e5, alpha=0.25, color='gold',
    label=r'D/H$\times 10^5 = 2.547 \pm 0.025$')
ax2.axhline(DH_obs*1e5, color='goldenrod', ls='--', lw=1)
ax2.axvspan(Neff_planck-2*Neff_planck_err, Neff_planck+2*Neff_planck_err,
    alpha=0.12, color='skyblue')
ax2.axvline(Neff_SM, color='green', lw=2, alpha=0.7)
ax2.axvline(Neff_s, color='blue', lw=2, ls='--', alpha=0.7)
ax2.axvline(Neff_ws, color='red', lw=2, ls=':', alpha=0.7)
ax2.set_xlabel(r'$N_{\rm eff}$')
ax2.set_ylabel(r'D/H $\times 10^5$')
ax2.set_xlim(2.5, 4.5); ax2.set_ylim(2.2, 3.2)
ax2.legend(loc='upper left', fontsize=9, framealpha=0.9)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, 'neff_bbn_three_scenarios.png'))
plt.close(fig)
print("  -> neff_bbn_three_scenarios.png")

print("\n" + "="*70)
print("PAPER 9 ANALYSIS - KEY NUMERICAL RESULTS")
print("="*70)
print(f"Delta_Neff (7 keV sterile) = {dN_sterile:.6f}")
print(f"Delta_Neff (1.08 TeV WIMP) = {dN_wimp:.8f}")
print(f"N_eff (SM + sterile) = {Neff_s:.6f}")
print(f"N_eff (SM + WIMP + sterile) = {Neff_ws:.6f}")
print("Higgsino: EXCLUDED (delta_LZ=377 < delta_lim=512-565 keV)")
print("Pseudo-Dirac: NOT excluded (Gamma_A/Gamma_lim ~ 5e-6)")
print("Neutron-philic O4: NOT excluded (C_sun/2Gamma_lim ~ 0.003)")
