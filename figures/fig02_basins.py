"""Figure 2 of the paper: share of initial states that reach the target under the replicator
dynamics. Left: over (phi, b_R5), read from ../model/figH_basin.npy. Right: against phi for three
values of b_R5 (200 random initial states, seed 1; cached in data_figI.npy).

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import os
import numpy as np
from print_style import style, save, rr, TW, HERE, plt
import botppp_model as B

style()
P = dict(B.DEFAULT)
GOOD = (1, 1, 1)
basin = np.load(rr('figH_basin.npy'))

cache = os.path.join(HERE, 'data_figI.npy')
if os.path.exists(cache):
    curves = np.load(cache, allow_pickle=True).item()
else:
    rng = np.random.default_rng(1)
    nstart = 200
    starts_r = rng.uniform(0.05, 0.95, size=(nstart, 3))
    phis_f = np.linspace(0, 1, 51)
    curves = {'phi': phis_f}
    for br in (0.0, 1.1, 1.9):
        X0 = np.tile(starts_r, (len(phis_f), 1))
        qv = {**P, 'bR5': br, 'phi': np.repeat(phis_f, nstart)}
        S = B.rk4_batch(X0, qv, T=2500, dt=0.1)
        good = ((np.abs(S - np.array(GOOD)) < 0.05).all(axis=1)).reshape(len(phis_f), nstart)
        curves[br] = good.mean(axis=1)
    np.save(cache, curves, allow_pickle=True)

fig = plt.figure(figsize=(TW, 2.3))
axL = fig.add_axes([0.095, 0.19, 0.30, 0.76])
cax = fig.add_axes([0.405, 0.19, 0.015, 0.76])
axR = fig.add_axes([0.625, 0.19, 0.355, 0.76])

im = axL.imshow(basin.T, origin='lower', extent=[0, 1, 0, 2], aspect='auto',
                cmap='viridis', vmin=0, vmax=1)
ph = np.linspace(P['cR'] / P['bR'] + 1e-6, 1, 200)
axL.plot(ph, P['bR'] - P['cR'] / ph, 'w-', lw=1.8)
axL.text(0.30, 1.52, '$\\phi^{*}$', color='w', fontsize=9.5, ha='center', va='center')
axL.set_xlabel('$\\phi$')
axL.set_ylabel('$b_{R5}$')
axL.set_ylim(0, 2)
cb = fig.colorbar(im, cax=cax)
cb.set_label('share reaching the target')
cb.set_ticks([0, 0.5, 1])

phis_f = curves['phi']
for br, col in [(0.0, 'tab:blue'), (1.1, 'tab:orange'), (1.9, 'tab:red')]:
    pstar = P['cR'] / (P['bR'] - br)
    lab = f'$b_{{R5}}={br:g}$'
    axR.plot(phis_f, curves[br], color=col, lw=1.6, label=lab)
    if pstar <= 1:
        axR.axvline(pstar, color=col, ls=':', lw=1.2, alpha=0.8)
axR.set_xlabel('$\\phi$')
axR.set_ylabel('share reaching the target')
axR.set_xlim(0, 1)
axR.set_ylim(-0.04, 1.04)
axR.grid(alpha=0.25)
axR.legend(loc='center right', bbox_to_anchor=(1.0, 0.70))
save(fig, 'fig2_basins.png')
