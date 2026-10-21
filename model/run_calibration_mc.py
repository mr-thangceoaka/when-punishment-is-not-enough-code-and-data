"""SI figure figS_calib_mc.png: Monte Carlo propagation of the calibration brackets.

Each quantity is drawn independently and uniformly from its bracket (lambda fixed, v ~ U(0, 0.1)).
Prints the share of draws in which the corridor lies above phi* and below phi_s and phi**.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import calibration_2024 as SL

CR_PRE, CR_POST = SL.CR_PRE, SL.CR_POST
LAM = round(SL.LAM, 3)
R_TARIFF = (round(SL.R_TARIFF[0], 3), round(SL.R_TARIFF[1], 3))
PHI_EMP = (round(SL.PHI_EMP[0], 2), round(SL.PHI_EMP[1], 2))
N = 200_000
rng = np.random.default_rng(1)
crpo, crpr = rng.uniform(*CR_POST, N), rng.uniform(*CR_PRE, N)
v = rng.uniform(0, 0.1, N)
D = 1 - rng.uniform(*R_TARIFF, N)
phi = rng.uniform(*PHI_EMP, N)
ps_pre, ps_post = crpr / D, crpo / D
phs = (2 * crpo + v + LAM) / (2 * D)
pss = (crpo + v + LAM) / D

p_pre = np.mean(phi > ps_pre)
p_post = np.mean(phi > ps_post)
p_band = np.mean((phi > ps_post) & (phi < phs))
print(f'N = {N}')
print(f'P[target stable, audit cost before ETC] = {p_pre:.4f}')
print(f'P[target stable, after ETC]             = {p_post:.4f}')
print(f'P[phi* < phi_emp < phi_s, after ETC]    = {p_band:.4f}')
print(f'P[phi_emp < phi**]                      = {np.mean(phi < pss):.4f}')
print(f'phi_s range in draws: {phs.min():.3f} - {phs.max():.3f};  phi** range: {pss.min():.3f} - {pss.max():.3f}')

FS = 14
plt.rcParams.update({'font.size': FS, 'axes.labelsize': FS + 1, 'xtick.labelsize': FS - 1, 'mathtext.fontset': 'dejavusans'})
fig, ax = plt.subplots(figsize=(7.0, 3.5))
bins = np.linspace(0, 1.6, 161)
ax.axvspan(*PHI_EMP, color='#cdebc9', label='proxy $\\phi_{\\rm emp}$ (bracket)')
ax.hist(ps_pre, bins, color='#a9c8e6', label='$\\phi^*$ before ETC')
ax.hist(ps_post, bins, color='#2f6fb0', label='$\\phi^*$ after ETC')
ax.hist(phs, bins, color='#9c9080', label='$\\phi_s$')
ax.hist(pss, bins, color='#e8877c', label='$\\phi^{**}$')
ax.axvline(1, color='0.4', lw=1)
ax.set_xlim(0, 1.6)
ax.set_yticks([])
ax.set_ylabel('draws')
ax.set_xlabel('$\\phi$')
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.24), ncol=3, frameon=False, fontsize=FS - 1,
          handlelength=1.2, columnspacing=1.2)
fig.savefig('figS_calib_mc.png', dpi=220, bbox_inches='tight')
plt.close(fig)
print('saved figS_calib_mc.png')
