"""Extended precision recomputation of the finite population checks of run_review_checks.py
(rule shape, reward size, retained fines, unequal population sizes) and of the beta = 10 panel.
Writes mp_checks.txt and a draft figure figM_adaptive.png.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import botppp_model as B
import fixation_general as G
from run_grid_mp import rho_mp, STATES, IDX, gth

mp.mp.dps = 60
P = dict(B.DEFAULT, bR5=0.0)

def pi_all(p, Zs=(50, 50, 50), beta=0.5, rule='punish', y0=0.5, s=0.05, kappa=0.0):
    L = mp.zeros(8, 8)
    q = {**B.DEFAULT, **p}
    for st in STATES:
        i = IDX[st]; x, y, z = st
        for pop in range(3):
            t = list(st); t[pop] = 1 - t[pop]; t = tuple(t)
            Z = Zs[pop]
            if pop == 1:
                mut = t[1]; dfs = []
                for j in range(1, Z):
                    yy = j/Z if mut == 1 else 1 - j/Z
                    dfs.append(G.conc_delta(p, x, z, yy, mut, G.weight(rule, yy, y0, s)))
            else:
                w = G.weight(rule, float(y), y0, s)
                T = dict(B.payoff_cells(p, w=w))
                if kappa:
                    a, c, r = T[(1, 0, 1)]; T[(1, 0, 1)] = (a, c, r + kappa*w*(q['um'] + q['uh']))
                    a, c, r = T[(0, 0, 1)]; T[(0, 0, 1)] = (a, c, r + kappa*w*q['um'])
                dfs = [T[t][pop] - T[st][pop]]*(Z - 1)
            L[i, IDX[t]] = rho_mp(dfs, beta)/3
    return gth(L)

def pt(p, **kw): return float(pi_all(p, **kw)[IDX[(1, 1, 1)]])

def lift(phis, **kw):
    for ph in phis:
        if pt(dict(P, phi=float(ph)), **{k: v for k, v in kw.items() if k != 'p'}) > 0.5:
            return float(ph)
    return np.nan

out = []
def log(t): print(t, flush=True); out.append(t)
phis = np.round(np.linspace(0, 1, 101), 2)

log('=== rule shape (logistic), Z=50 ===')
for beta in (0.5, 10.0):
    for y0 in (0.1, 0.25, 0.5, 0.75, 0.9):
        row = [lift(phis, beta=beta, rule='logistic', y0=y0, s=s) for s in (0.01, 0.05, 0.1, 0.2)]
        log(f'beta={beta:4.1f} y0={y0:4.2f}: ' + '  '.join(f'{v:.2f}' for v in row))

log('=== reward size r = k u, step rule y0=1/2, Z=50 ===')
for k in (0.4, 0.6, 0.7, 0.8, 1.0, 1.25, 1.5):
    def liftk(beta):
        for ph in phis:
            if pt(dict(P, phi=float(ph), rm=0.8*k, rh=0.8*k), beta=beta, rule='step') > 0.5:
                return float(ph)
        return np.nan
    log(f'k={k:4.2f}: beta=0.5 {liftk(0.5):.2f}  beta=10 {liftk(10.0):.2f}')

log('=== retained fines kappa, punishment, Z=50 ===')
for kappa in (0.0, 0.1, 0.25, 0.5, 1.0):
    log(f'kappa={kappa:4.2f}: beta=0.5 {lift(phis, beta=0.5, rule="punish", kappa=kappa):.2f}  '
        f'beta=10 {lift(phis, beta=10.0, rule="punish", kappa=kappa):.2f}')

log('=== unequal population sizes (Z_U, Z_C, Z_S) ===')
for Zs in [(50, 50, 50), (1000, 50, 50), (1000, 10, 20), (1000, 100, 10), (200, 5, 5), (10, 100, 100)]:
    r = [lift(phis, Zs=Zs, beta=b, rule=rule) for b in (0.5, 10.0) for rule in ('punish', 'step')]
    log(f'{Zs}: beta=0.5 punish {r[0]:.2f} rule {r[1]:.2f} | beta=10 punish {r[2]:.2f} rule {r[3]:.2f}')

log('=== figM curves ===')
curves = {}
for beta in (0.5, 10.0):
    for rule in ('punish', 'reward', 'step'):
        curves[(rule, beta)] = np.array([pt(dict(P, phi=float(ph)), beta=beta, rule=rule) for ph in phis])
        pg = curves[(rule, beta)]
        log(f'beta={beta} {rule}: lift off {phis[np.argmax(pg > 0.5)] if pg.max() > 0.5 else np.nan:.2f}')
open('mp_checks.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')

plt.rcParams.update({'font.size': 13, 'axes.titlesize': 15, 'axes.labelsize': 16,
                     'xtick.labelsize': 12.5, 'ytick.labelsize': 12.5, 'lines.linewidth': 2.4})
NAME = {'punish': 'punishment', 'reward': 'reward', 'step': 'state dependent rule'}
COL = {'punish': 'tab:red', 'reward': 'tab:orange', 'step': 'tab:green'}
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
for ax, beta in zip(axes, (0.5, 10.0)):
    for rule in ('punish', 'reward', 'step'):
        ax.plot(phis, curves[(rule, beta)], color=COL[rule], label=NAME[rule])
    ax.axvline(0.10, color='gray', ls=':', lw=2); ax.axvline(0.50, color='gray', ls='--', lw=2)
    ax.set_xlabel('$\\phi$'); ax.set_title(f'Z=50, $\\beta$={beta}')
axes[0].set_ylabel('$\\pi$(target state)')
fig.tight_layout()
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=3, fontsize=13, framealpha=0.95, bbox_to_anchor=(0.5, -0.06))
fig.savefig('figM_adaptive.png', dpi=200, bbox_inches='tight'); plt.close(fig)
print('wrote mp_checks.txt figM_adaptive.png')
