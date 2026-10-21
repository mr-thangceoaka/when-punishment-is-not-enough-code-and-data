"""SI figures at print size: figSI_plane.png (vector field on the face x = 1), figSI_grid.png
(switching point over population size and selection strength, from ../model/figN_tab_*.npy written
by run_grid_mp.py) and figSI_rule.png (robustness of the state dependent rule; numbers read from
../model/review_checks.txt).

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import re
import numpy as np
from matplotlib.lines import Line2D
from print_style import style, save, rr, TW_SI, plt
import botppp_model as B

style()
P = dict(B.DEFAULT)



def _field(Y, Z, q):
    _, dC, dR = B.deltas_vec(1.0, Y, Z, q)
    return Y * (1 - Y) * dC, Z * (1 - Z) * dR


def _kind(y, z, q, h=1e-5):
    def F(a, b):
        _, dC, dR = B.deltas_vec(1.0, a, b, q)
        return np.array([a * (1 - a) * dC, b * (1 - b) * dR])
    J = np.column_stack([(F(y + h, z) - F(y - h, z)) / (2 * h), (F(y, z + h) - F(y, z - h)) / (2 * h)])
    ev = np.linalg.eigvals(J).real
    return 'stable' if np.all(ev < -1e-9) else 'unstable' if np.all(ev > 1e-9) else 'saddle'


g = np.linspace(0, 1, 400)
Yg, Zg = np.meshgrid(g, g)
fig = plt.figure(figsize=(0.95 * TW_SI, 3.25))
axes = [fig.add_axes([0.075, 0.25, 0.345, 0.66]), fig.add_axes([0.53, 0.25, 0.345, 0.66])]
cax = fig.add_axes([0.895, 0.25, 0.015, 0.66])
strm = None
for ax, (lab, q) in zip(axes, [(r'valve open ($b_{R5}=0$)', dict(P, phi=1.0, bR5=0.0)),
                               (r'valve closed ($b_{R5}=1.9$)', dict(P, phi=1.0, bR5=1.9))]):
    U, V = _field(Yg, Zg, q)
    spd = np.hypot(U, V)
    strm = ax.streamplot(g, g, U, V, color=spd, cmap='viridis', density=0.8, linewidth=0.7,
                         arrowsize=0.7, broken_streamlines=False, zorder=1)
    zc = (q['cP'] + q['gm']) / (q['um'] + q['uh'] - q['gh'])
    if 0 < zc < 1:
        ax.axhline(zc, ls=(0, (5, 3)), color='crimson', lw=1.1, zorder=5)
    kw = dict(zorder=7, clip_on=False)
    for fy, fz in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        k = _kind(fy, fz, q)
        ax.scatter([fy], [fz], s=70, marker='o' if k != 'saddle' else 's', facecolor='white', edgecolor='none', **kw)
        if k == 'stable':
            halo = 'seagreen' if (fy, fz) == (1, 1) else 'firebrick'
            ax.scatter([fy], [fz], s=190, facecolor=halo, alpha=.28, edgecolor='none', zorder=6, clip_on=False)
            ax.scatter([fy], [fz], s=55, marker='o', facecolor='black', edgecolor='white', linewidths=1.0, **kw)
        elif k == 'unstable':
            ax.scatter([fy], [fz], s=50, marker='o', facecolor='white', edgecolor='black', linewidths=1.0, **kw)
        else:
            ax.scatter([fy], [fz], s=42, marker='s', facecolor='0.82', edgecolor='black', linewidths=0.9, **kw)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r'$y$ (fraction Comply)')
    ax.set_ylabel(r'$z$ (fraction Audit)')
    ax.set_title(lab, pad=4)
cb = fig.colorbar(strm.lines, cax=cax)
cb.set_label(r'flow speed $\|(\dot y,\dot z)\|$')
fig.legend(handles=[
    Line2D([], [], marker='o', color='w', mfc='black', mec='white', ms=6, label='stable (attractor)'),
    Line2D([], [], marker='o', color='w', mfc='white', mec='black', ms=6, label='unstable'),
    Line2D([], [], marker='s', color='w', mfc='0.82', mec='black', ms=6, label='saddle'),
    Line2D([], [], ls=(0, (5, 3)), color='crimson', label=r'$\dot y=0$ nullcline')],
    loc='lower center', bbox_to_anchor=(0.47, 0.0), ncol=4, frameon=False)
save(fig, 'figSI_plane.png')

Zs_list = [25, 50, 100, 200, 500]
betas = [0.1, 0.5, 2.0, 10.0]
tab = {'punish': np.load(rr('figN_tab_punish.npy')), 'step': np.load(rr('figN_tab_adaptive.npy'))}


def heat_pair(tabs, titles, xt, yt, xlabel, ylabel, vmax, name):
    fig = plt.figure(figsize=(TW_SI, 2.55))
    axs = [fig.add_axes([0.075, 0.17, 0.37, 0.72]), fig.add_axes([0.52, 0.17, 0.37, 0.72])]
    cax = fig.add_axes([0.905, 0.17, 0.015, 0.72])
    for ax, M, ttl in zip(axs, tabs, titles):
        im = ax.imshow(M, cmap='viridis', vmin=0, vmax=vmax, aspect='auto')
        ax.set_xticks(range(len(xt)))
        ax.set_xticklabels(xt)
        ax.set_yticks(range(len(yt)))
        ax.set_yticklabels(yt)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(ttl, pad=4)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                ax.text(j, i, f'{M[i, j]:.2f}', ha='center', va='center', color='w', fontsize=8,
                        fontweight='bold')
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('switching point ($\\phi$)')
    save(fig, name)


heat_pair([tab['punish'], tab['step']],
          ['punishment: tends to 0.30 as $\\beta Z$ grows', 'state dependent rule: tends to 0.1'],
          betas, Zs_list, '$\\beta$', '$Z$', 0.7, 'figSI_grid.png')

txt = open(rr('review_checks.txt'), encoding='utf-8').read()
blk = txt[txt.index('=== R1b'):txt.index('=== R1c')]
rows = re.findall(r'y0=([0-9.]+):\s+s=0\.01:([0-9.]+)\s+s=0\.05:([0-9.]+)\s+s=0\.1:([0-9.]+)\s+s=0\.2:([0-9.]+)', blk)
assert len(rows) == 10, rows
A = np.array([[float(v) for v in r[1:]] for r in rows])
Y0 = [float(r[0]) for r in rows[:5]]
S = [0.01, 0.05, 0.1, 0.2]
heat_pair([A[:5], A[5:]], ['$Z=50$, $\\beta=0.5$', '$Z=50$, $\\beta=10$'], S, Y0,
          'steepness $s$', 'switching point $y_0$', 0.6, 'figSI_rule.png')
