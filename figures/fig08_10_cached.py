"""Figure 3 of the paper (fig8_rule.png) and SI figures fig9_users.png and fig10_welfare.png,
drawn from cached results in ../model: figMQ_curves.npy, basins.npy, ext_E1.npy and ext_E4.npy.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
from print_style import style, save, rr, TW, plt
import botppp_model as B

style()

from matplotlib.lines import Line2D
phis = np.round(np.linspace(0, 1, 101), 2)
curves = np.load(rr('figMQ_curves.npy'), allow_pickle=True).item()
d = np.load(rr('basins.npy'), allow_pickle=True).item()
COL = {'punish': 'tab:red', 'reward': 'tab:orange', 'step': 'tab:green'}
fig = plt.figure(figsize=(TW, 2.6))
W, H, Y0 = 0.235, 0.60, 0.165
axes = [fig.add_axes([0.085, Y0, W, H]), fig.add_axes([0.36, Y0, W, H]), fig.add_axes([0.745, Y0, W, H])]
for ax, beta, ttl in zip(axes[:2], (0.5, 10.0), (r'(a) finite, $\beta=0.5$', r'(b) finite, $\beta=10$')):
    for rule in ('punish', 'reward', 'step'):
        ax.plot(phis, curves[(rule, beta)], color=COL[rule], lw=1.5)
    ax.axvline(0.1, color='gray', ls=':', lw=1.0)
    ax.axvline(0.5, color='gray', ls='--', lw=1.0)
    ax.set_xlabel(r'$\phi$')
    ax.set_title(ttl, pad=3)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.03, 1.03)
axes[0].set_ylabel('long run share of target')
axes[1].set_yticklabels([])
ax = axes[2]
sty = {'punish': ('tab:red', 'o'), 'y0=0.5': ('tab:green', 's'), 'y0=0.333': ('tab:blue', 'D'),
       'y0=0.2': ('tab:purple', '^')}
for k in [k for k in d if k != 'phi']:
    c, mk = sty[k]
    ax.plot(d['phi'], d[k], marker=mk, ms=3, color=c, lw=1.2)
ax.axvline(0.1, color='gray', ls=':', lw=1.0)
ax.axvline(0.5, color='gray', ls='--', lw=1.0)
ax.axvline(0.367, color='tab:blue', ls='-.', lw=1.0)
ax.set_xlabel(r'$\phi$')
ax.set_ylabel('share of starts\nreaching target')
ax.set_title('(c) replicator', pad=3)
ax.set_xlim(0, 0.8)
ax.set_ylim(-0.03, 1.03)
leg = [Line2D([0], [0], color='tab:red', lw=1.5, label='punishment'),
       Line2D([0], [0], color='tab:orange', lw=1.5, label='reward'),
       Line2D([0], [0], color='tab:green', lw=1.5, label='state dependent rule, $y_0=1/2$'),
       Line2D([0], [0], color='tab:blue', lw=1.5, label='rule, $y_0=1/3$ (c only)'),
       Line2D([0], [0], color='tab:purple', lw=1.5, label='rule, $y_0=0.2$ (c only)')]
fig.legend(handles=leg, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False,
           columnspacing=1.0)
save(fig, 'fig8_rule.png')

E1 = np.load(rr('ext_E1.npy'), allow_pickle=True).item()
PO = np.round(np.arange(0.0, 0.81, 0.05), 2)
PQ = np.round(np.arange(0.0, 0.81, 0.05), 2)
P = dict(B.DEFAULT, bR5=0.0)
q_c = lambda po: (P['cP'] + P['gm'] + P['gh'] * (1 - po)) / (P['bP'] + P['gh'])
fig = plt.figure(figsize=(TW, 2.2))
W, H, Y0 = 0.245, 0.58, 0.20
axs = [fig.add_axes([0.075, Y0, W, H]), fig.add_axes([0.36, Y0, W, H]), fig.add_axes([0.645, Y0, W, H])]
cax = fig.add_axes([0.905, Y0, 0.014, H])
fig.text((0.075 + 0.36 + W) / 2, 0.955, 'under punishment', ha='center', va='center', fontsize=9)
fig.text(0.645 + W / 2, 0.955, 'under the rule', ha='center', va='center', fontsize=9)
for xa, xb in [(0.075, 0.36 + W), (0.645, 0.645 + W)]:
    fig.add_artist(plt.Line2D([xa, xb], [0.915, 0.915], color='0.4', lw=0.7))
for k, (ax, (A, ttl)) in enumerate(zip(axs, [(E1['punish'][:, :, 0], 'target state'),
                                             (E1['punish'][:, :, 1], '(Use, Comply, Lax)'),
                                             (E1['step'][:, :, 0] + E1['step'][:, :, 1],
                                              'compliance')])):
    im = ax.imshow(A, origin='lower', cmap='viridis', vmin=0, vmax=1, aspect='auto',
                   extent=[PQ[0] - 0.025, PQ[-1] + 0.025, PO[0] - 0.025, PO[-1] + 0.025])
    ax.plot([q_c(po) for po in PO], PO, 'w--', lw=1.2)
    ax.axhline(0.1, color='w', ls=':', lw=1.1)
    ax.set_xlabel('$\\phi_q$')
    if k == 0:
        ax.set_ylabel('$\\phi$')
    else:
        ax.set_yticklabels([])
    ax.set_title(ttl, pad=3)
    ax.set_xticks([0, 0.4, 0.8])
    ax.set_yticks([0, 0.4, 0.8])
cb = fig.colorbar(im, cax=cax)
cb.set_label('long run share')
cb.set_ticks([0, 0.5, 1])
save(fig, 'fig9_users.png')

E4 = np.load(rr('ext_E4.npy'), allow_pickle=True).item()
PH4 = E4['phi']
lab = {'none': 'no audit', 'punish': 'punishment', 'logistic': 'state dependent rule ($y_0=1/3$)'}
col = {'none': 'gray', 'punish': 'tab:red', 'logistic': 'tab:blue'}
fig = plt.figure(figsize=(TW, 2.25))
a0 = fig.add_axes([0.10, 0.19, 0.37, 0.64])
a1 = fig.add_axes([0.61, 0.19, 0.37, 0.64])
for r in ('none', 'punish', 'logistic'):
    a0.plot(PH4, E4[r][:, 0], marker='o', ms=3, color=col[r], lw=1.4, label=lab[r])
    a1.plot(PH4, E4[r][:, 1], marker='o', ms=3, color=col[r], lw=1.4, label=lab[r])
for ax in (a0, a1):
    ax.axvline(0.1, color='gray', ls=':', lw=1.0)
    ax.axvline(0.5, color='gray', ls='--', lw=1.0)
    ax.set_xlabel('$\\phi$')
a0.set_ylabel('discounted welfare per period')
a1.set_ylabel('discounted cost of\nrewards per period')
h, l = a0.get_legend_handles_labels()
fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False, columnspacing=1.2)
save(fig, 'fig10_welfare.png')
