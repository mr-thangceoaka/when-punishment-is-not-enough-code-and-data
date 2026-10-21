"""Figure 4 of the paper: the case-study corridor in the plane of the toll ratio b_R5/b_R and phi.

Units: b_R = 1 (expressway revenue per vehicle in 2024). Brackets from calibration_2024.py:
  c_R:         0.06-0.12 before electronic toll collection (ETC), 0.02-0.045 after (1 June 2022)
  v + lambda:  0.593-0.693 (lambda = 2024 interest / 2024 expressway revenue; 0 <= v <= 0.1)
  b_R5/b_R:    0.407-0.462 over the five vehicle classes (full-trip fares from 29 December 2023)
  phi_emp:     0.25-0.33 (operator's 2016 count)
Writes fig6_calibration_data.png and prints the threshold brackets.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import calibration_2024 as SL
CR_PRE = SL.CR_PRE
CR_POST = SL.CR_POST
VLAM = (round(SL.LAM, 3), round(SL.LAM, 3) + 0.1)
R_TARIFF = (round(SL.R_TARIFF[0], 3), round(SL.R_TARIFF[1], 3))
R_REV = round(SL.R_REALISED, 3)
PHI_EMP = (round(SL.PHI_EMP[0], 2), round(SL.PHI_EMP[1], 2))
FS = 15
plt.rcParams.update({'font.size': FS, 'axes.labelsize': FS + 1, 'xtick.labelsize': FS - 1,
                     'ytick.labelsize': FS - 1, 'mathtext.fontset': 'dejavusans'})

r = np.linspace(0, 0.75, 600)
D = 1 - r
mid = lambda a: 0.5 * (a[0] + a[1])
TOP = 1.0
clip = lambda y: np.minimum(y, TOP)

ps_lo, ps_hi = CR_POST[0] / D, CR_POST[1] / D
pre_lo, pre_hi = CR_PRE[0] / D, CR_PRE[1] / D
pss_lo, pss_hi = (CR_POST[0] + VLAM[0]) / D, (CR_POST[1] + VLAM[1]) / D
phs_lo = (2 * CR_POST[0] + VLAM[0]) / (2 * D)
phs_hi = (2 * CR_POST[1] + VLAM[1]) / (2 * D)
ps_mid, pss_mid = mid(CR_POST) / D, (mid(CR_POST) + mid(VLAM)) / D

fig = plt.figure(figsize=(7.2, 5.3))
ax = fig.add_axes([0.12, 0.13, 0.60, 0.76])
C_UNST, C_TRAP, C_REACH = '#dedede', '#fde9d2', '#cdebc9'
ax.fill_between(r, 0, clip(ps_mid), color=C_UNST, lw=0)
ax.fill_between(r, clip(ps_mid), clip(pss_mid), color=C_TRAP, lw=0)
ax.fill_between(r, clip(pss_mid), TOP, where=pss_mid < TOP, color=C_REACH, lw=0)
ax.fill_between(r, ps_lo, ps_hi, color='#6f9fcf', lw=0)
ax.plot(r, clip(pre_lo), color='tab:blue', ls='--', lw=1.8)
ax.plot(r, clip(pre_hi), color='tab:blue', ls='--', lw=1.8)
ax.fill_between(r, clip(phs_lo), clip(phs_hi), color='#bfb5a8', lw=0)
ax.fill_between(r, clip(pss_lo), clip(pss_hi), color='#e8877c', lw=0)
ax.add_patch(Rectangle((R_TARIFF[0], PHI_EMP[0]), R_TARIFF[1] - R_TARIFF[0], PHI_EMP[1] - PHI_EMP[0],
                       fill=True, fc=(0.0, 0.5, 0.0, 0.18), ec='darkgreen', lw=2.6, zorder=5))
ax.annotate('Hanoi Hai Phong\ncorridor', xy=(R_TARIFF[0], mid(PHI_EMP)), xytext=(0.035, 0.235),
            fontsize=FS, color='darkgreen', va='center',
            arrowprops=dict(arrowstyle='-', color='darkgreen', lw=1.4))

ax.text(0.585, 0.64, 'trapped\nband', fontsize=FS + 2, ha='center', va='center', fontweight='bold', color='#9a5b13')
ax.text(0.050, 0.555, '$\\phi_s$', fontsize=FS + 1, color='0.25')
ax.text(0.215, 0.70, '$\\phi^{**}$', fontsize=FS + 1, color='firebrick')
ax.set_xlim(0, 0.75)
ax.set_ylim(0, TOP)
ax.set_xlabel('toll ratio $b_{R5}/b_R$')
ax.set_ylabel('discerning fraction of users, $\\phi$')

xr = 1.02
def right_label(y, text, color):
    ax.text(xr, y, text, transform=ax.get_yaxis_transform(), fontsize=FS, color=color, va='center', clip_on=False)
right_label(clip(pre_hi[-1] + pre_lo[-1]) / 2 + 0.03, '$\\phi^*$ before ETC', 'tab:blue')
right_label((ps_hi[-1] + ps_lo[-1]) / 2, '$\\phi^*$ after ETC', 'tab:blue')
right_label(0.035, 'target not stable', '0.30')
ax.annotate('', xy=(1.07, (ps_hi[-1] + ps_lo[-1]) / 2 + 0.05), xycoords=ax.get_yaxis_transform(),
            xytext=(1.07, (pre_hi[-1] + pre_lo[-1]) / 2 - 0.01), textcoords=ax.get_yaxis_transform(),
            arrowprops=dict(arrowstyle='->', color='tab:blue', lw=1.8), annotation_clip=False)
ax.text(1.10, ((ps_hi[-1] + ps_lo[-1]) / 2 + (pre_hi[-1] + pre_lo[-1]) / 2) / 2 + 0.02, 'ETC 2022',
        transform=ax.get_yaxis_transform(), fontsize=FS - 1, color='tab:blue', va='center', clip_on=False)
ax.text(0.0, 1.035, 'reachable from every state ($\\phi > \\phi^{**}$)', transform=ax.transAxes,
        fontsize=FS, color='darkgreen', va='bottom', clip_on=False)
ax.annotate('', xy=(0.03, 0.965), xytext=(0.06, 1.03), textcoords=ax.transAxes, xycoords='data',
            arrowprops=dict(arrowstyle='->', color='darkgreen', lw=1.4), annotation_clip=False)

fig.savefig('fig6_calibration_data.png', dpi=220)
plt.close(fig)

one_r = (1 - R_TARIFF[1], 1 - R_TARIFF[0])
print(f"1-r (tariff)          : {one_r[0]:.3f} - {one_r[1]:.3f}; revenue basis {1-R_REV:.3f}")
print(f"phi*  before ETC      : {CR_PRE[0]/one_r[1]:.3f} - {CR_PRE[1]/one_r[0]:.3f}")
print(f"phi*  after ETC       : {CR_POST[0]/one_r[1]:.3f} - {CR_POST[1]/one_r[0]:.3f}")
print(f"phi** punish after ETC: {(CR_POST[0]+VLAM[0])/one_r[1]:.3f} - {(CR_POST[1]+VLAM[1])/one_r[0]:.3f}")
print(f"phi_s                 : {(2*CR_POST[0]+VLAM[0])/(2*one_r[1]):.3f} - {(2*CR_POST[1]+VLAM[1])/(2*one_r[0]):.3f}")
print(f"phi_emp               : {PHI_EMP[0]:.2f} - {PHI_EMP[1]:.2f}")
