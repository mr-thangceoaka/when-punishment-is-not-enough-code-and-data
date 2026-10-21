"""Switching point over (Z, beta) in extended precision (mpmath, 60 digits).

At large beta*Z the fixation probabilities fall far below the double precision range, and the
stationary distribution of the 8 x 8 transition matrix becomes unreliable in double precision.
Here the fixation probability (general formula, rule evaluated at every mutant count) and the
stationary distribution (algorithm of Grassmann, Taksar and Heyman) are computed in extended
precision. Writes figN_tab_punish.npy, figN_tab_adaptive.npy, grid_mp.txt and a draft figure
figN_robustness.png.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import sys
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import botppp_model as B
import fixation_general as G

mp.mp.dps = 60
P = dict(B.DEFAULT, bR5=0.0)
STATES = B.STATES
IDX = {s: i for i, s in enumerate(STATES)}

def rho_mp(dfs, beta):
    acc = mp.mpf(0); tot = mp.mpf(1)
    for d in dfs:
        acc += -mp.mpf(beta)*mp.mpf(d)
        tot += mp.e**acc
    return 1/tot


def gth(Lo):
    """Stationary distribution by the Grassmann, Taksar and Heyman algorithm.
    Uses only off diagonal rates, so no cancellation occurs however small they are."""
    n = Lo.rows
    A = Lo.copy()
    for k in range(n - 1, 0, -1):
        s = sum(A[k, j] for j in range(k))
        for i in range(k):
            A[i, k] = A[i, k]/s
        for i in range(k):
            for j in range(k):
                if i != j:
                    A[i, j] += A[i, k]*A[k, j]
    pi = [mp.mpf(0)]*n
    pi[0] = mp.mpf(1)
    for k in range(1, n):
        pi[k] = sum(pi[i]*A[i, k] for i in range(k))
    tot = sum(pi)
    return [v/tot for v in pi]

def pi_target(p, Zs, beta, rule):
    L = mp.zeros(8, 8)
    for st in STATES:
        i = IDX[st]; x, y, z = st
        for pop in range(3):
            t = list(st); t[pop] = 1 - t[pop]; t = tuple(t)
            Z = Zs[pop]
            if pop == 1:
                mut = t[1]; dfs = []
                for j in range(1, Z):
                    yy = j/Z if mut == 1 else 1 - j/Z
                    dfs.append(G.conc_delta(p, x, z, yy, mut, G.weight(rule, yy)))
            else:
                T = B.payoff_cells(p, w=G.weight(rule, float(y)))
                dfs = [T[t][pop] - T[st][pop]]*(Z - 1)
            L[i, IDX[t]] = rho_mp(dfs, beta)/3
    pi = gth(L)
    return float(pi[IDX[(1, 1, 1)]])

def liftoff(Zs, beta, rule, phis):
    for ph in phis:
        if pi_target(dict(P, phi=float(ph)), Zs, beta, rule) > 0.5:
            return float(ph)
    return np.nan

if __name__ == '__main__':
    phis = np.round(np.linspace(0, 1, 51), 2)
    Zs_list = [25, 50, 100, 200, 500]
    betas = [0.1, 0.5, 2.0, 10.0]
    out = []
    tab = {}
    for rule in ('punish', 'step'):
        M = np.zeros((len(Zs_list), len(betas)))
        for i, Z in enumerate(Zs_list):
            for j, b in enumerate(betas):
                M[i, j] = liftoff((Z, Z, Z), b, rule, phis)
            line = f'{rule:7s} Z={Z:4d}: ' + '  '.join(f'{M[i,j]:.2f}' for j in range(len(betas)))
            print(line, flush=True); out.append(line)
        tab[rule] = M
    line = f'check Z=50 beta=0.5: punish {liftoff((50,50,50),0.5,"punish",phis):.2f} step {liftoff((50,50,50),0.5,"step",phis):.2f}'
    print(line); out.append(line)
    open('grid_mp.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
    np.save('figN_tab_punish.npy', tab['punish']); np.save('figN_tab_adaptive.npy', tab['step'])

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.7), gridspec_kw={'wspace': 0.28})
    for ax, key, name, tgt in zip(axes, ('punish', 'step'), ('punishment', 'state dependent rule'), (0.5, 0.1)):
        M = tab[key]
        im = ax.imshow(M, cmap='viridis', vmin=0, vmax=0.7, aspect='auto')
        ax.set_xticks(range(len(betas))); ax.set_xticklabels(betas)
        ax.set_yticks(range(len(Zs_list))); ax.set_yticklabels(Zs_list)
        ax.set_xlabel('$\\beta$', fontsize=15); ax.set_ylabel('Z', fontsize=15)
        ax.set_title(f'{name}: tends to {tgt} as $\\beta Z$ grows', fontsize=13, pad=8)
        for i in range(len(Zs_list)):
            for j in range(len(betas)):
                ax.text(j, i, f'{M[i,j]:.2f}', ha='center', va='center', color='w', fontsize=12, fontweight='bold')
    cb = fig.colorbar(im, ax=axes, fraction=0.03); cb.set_label('lift off value of $\\phi$', fontsize=14)
    fig.suptitle('Lift off value of $\\phi$ over population size Z and selection strength $\\beta$',
                 fontsize=15.5, y=1.06)
    fig.savefig('figN_robustness.png', dpi=200, bbox_inches='tight'); plt.close(fig)
    print('wrote figN_robustness.png grid_mp.txt')
