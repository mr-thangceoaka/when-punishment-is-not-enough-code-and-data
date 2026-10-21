"""SI figure fig7_finite.png: long run share of the target in finite populations. Left: over
(phi, b_R5), read from ../model/figK_pigood.npy. Right: against phi (Z = 50, beta = 0.5).

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
from print_style import style, save, rr, TW, plt
import botppp_model as B

style()
P = dict(B.DEFAULT)
PG = np.load(rr('figK_pigood.npy'))

fig = plt.figure(figsize=(TW, 2.45))
axL = fig.add_axes([0.095, 0.18, 0.30, 0.70])
cax = fig.add_axes([0.405, 0.18, 0.015, 0.70])
axR = fig.add_axes([0.625, 0.18, 0.355, 0.70])

im = axL.imshow(PG.T, origin='lower', extent=[0, 1, 0, 2], aspect='auto', cmap='viridis', vmin=0, vmax=1)
ph = np.linspace(P['cR'] / P['bR'] + 1e-6, 1, 300)
axL.plot(ph, P['bR'] - P['cR'] / ph, c='cyan', lw=1.8)
k2 = P['cR'] + P['v'] + P['lam']
ph2 = np.linspace(k2 / P['bR'] + 1e-6, 1, 300)
axL.plot(ph2, P['bR'] - k2 / ph2, color='k', ls='--', lw=1.6)
axL.text(0.27, 1.55, '$\\phi^{*}$', color='cyan', fontsize=9.5, ha='center', va='center')
axL.text(0.80, 0.12, '$\\phi^{**}$', color='k', fontsize=9.5, ha='center', va='center')
axL.set_xlabel('$\\phi$')
axL.set_ylabel('$b_{R5}$')
axL.set_ylim(0, 2)
cb = fig.colorbar(im, cax=cax)
cb.set_label('$\\pi$(target state)')
cb.set_ticks([0, 0.5, 1])

phis3 = np.linspace(0, 1, 101)
for br, col in [(0.0, 'tab:blue'), (1.1, 'tab:orange')]:
    pg = [B.pi_good(dict(P, phi=ph_, bR5=br), Z=50, beta=0.5) for ph_ in phis3]
    ps = P['cR'] / (P['bR'] - br)
    pss = k2 / (P['bR'] - br)
    axR.plot(phis3, pg, color=col, lw=1.6, label=f'$b_{{R5}}={br:g}$')
    for val, ls in [(ps, ':'), (pss, '--')]:
        if val <= 1:
            axR.axvline(val, color=col, ls=ls, lw=1.2, alpha=0.8)
axR.set_xlabel('$\\phi$')
axR.set_ylabel('$\\pi$(target state)')
axR.set_xlim(0, 1)
axR.set_ylim(-0.04, 1.04)
axR.grid(alpha=0.25)
axR.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False)
save(fig, 'fig7_finite.png')
