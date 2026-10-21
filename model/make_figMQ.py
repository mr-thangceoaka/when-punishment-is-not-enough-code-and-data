"""Curves for the state dependent rule figure (Figure 3 of the paper), cached for the print figure.

(a), (b): long run share of the target state in finite populations for punishment, pure reward and
the step rule (extended precision, general fixation formula; Z = 50, beta = 0.5 and 10);
(c): share of initial states reaching the target under the replicator dynamics (from basins.npy).
Writes figMQ_curves.npy (read by figures/fig08_10_cached.py) and a draft figure figMQ.png.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import run_grid_mp as R

phis = np.round(np.linspace(0, 1, 101), 2)
cache = 'figMQ_curves.npy'
if os.path.exists(cache):
    curves = np.load(cache, allow_pickle=True).item()
else:
    curves = {}
    for beta in (0.5, 10.0):
        for rule in ('punish', 'reward', 'step'):
            curves[(rule, beta)] = np.array([R.pi_target(dict(R.P, phi=float(p)), (50, 50, 50), beta, rule)
                                             for p in phis])
    np.save(cache, curves, allow_pickle=True)
for k, v in curves.items():
    print(k, 'switching point', phis[np.argmax(v > 0.5)] if v.max() > 0.5 else None)

d = np.load('basins.npy', allow_pickle=True).item()
plt.rcParams.update({'font.size': 12, 'axes.labelsize': 14, 'axes.titlesize': 13})
fig, axes = plt.subplots(1, 3, figsize=(15, 3.9), gridspec_kw={'wspace': 0.28})
NAME = {'punish': 'punishment', 'reward': 'reward', 'step': 'state dependent rule'}
COL = {'punish': 'tab:red', 'reward': 'tab:orange', 'step': 'tab:green'}
for ax, beta, ttl in zip(axes[:2], (0.5, 10.0), ('(a) finite, $\\beta=0.5$', '(b) finite, $\\beta=10$')):
    for rule in ('punish', 'reward', 'step'):
        ax.plot(phis, curves[(rule, beta)], color=COL[rule], lw=2.4, label=NAME[rule])
    ax.axvline(0.1, color='gray', ls=':', lw=1.6); ax.axvline(0.5, color='gray', ls='--', lw=1.6)
    ax.set_xlabel('$\\phi$'); ax.set_title(ttl)
axes[0].set_ylabel('long run share of target')
axes[0].legend(fontsize=10.5, loc='center right', framealpha=0.95)
ax = axes[2]
sty = {'punish': ('tab:red', 'o', 'punishment'), 'y0=0.5': ('tab:green', 's', 'rule, $y_0=0.5$'),
       'y0=0.2': ('tab:purple', '^', 'rule, $y_0=0.2$')}
for k in [k for k in d if k != 'phi']:
    c, mk, lab = sty.get(k, ('tab:blue', 'D', 'rule, $y_0=1/3$'))
    ax.plot(d['phi'], d[k], marker=mk, ms=5, color=c, lw=2.0, label=lab)
ax.axvline(0.1, color='gray', ls=':', lw=1.6); ax.axvline(0.5, color='gray', ls='--', lw=1.6)
ax.axvline(0.367, color='tab:blue', ls='-.', lw=1.4)
ax.set_xlabel('$\\phi$'); ax.set_ylabel('share of starts reaching target')
ax.set_title('(c) replicator dynamics'); ax.set_ylim(-0.03, 1.03)
ax.legend(fontsize=10, loc='lower right', framealpha=0.95)
fig.savefig('figMQ.png', dpi=200, bbox_inches='tight'); plt.close(fig)
print('saved figMQ.png')
