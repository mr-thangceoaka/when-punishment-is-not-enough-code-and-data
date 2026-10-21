"""SI figure fig6_cube.png: stationary distributions of the finite population Markov chain
(Z = 50, beta = 0.5) for four mechanisms, drawn on the cube of the eight pure states.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
from matplotlib.patches import FancyArrowPatch, Circle
from print_style import style, save, TW, plt
import botppp_model as B

style()
P = dict(bU=1.0, bP=1.0, cP=0.3, gm=0.4, gh=0.4, um=0.8, uh=0.8,
         bR=2.0, bR5=0.0, cR=0.2, v=0.4, lam=0.4, theta=0.3)
Z, BETA = 50, 0.5
C_GOOD = ('#8ecf8e', '#2e7d32', '#1b4d1e')
C_BAD = ('#e8938f', '#c0392b', '#7d221b')
C_MID = ('#c3d9ee', '#5a86b8', '#2f4a68')
A_STR, A_WEAK = '#c0392b', '#8a949d'


def chain(q):
    idx = {s: i for i, s in enumerate(B.STATES)}
    L = np.zeros((8, 8))
    for s in B.STATES:
        i = idx[s]
        f = B.fitness(q, [float(v) for v in s], sharp=True)
        for pop in range(3):
            t = list(s); t[pop] = 1 - t[pop]; t = tuple(t)
            L[i, idx[t]] = B.fixation(f[pop][t[pop]] - f[pop][s[pop]], Z, BETA) / 3.0
        L[i, i] = 1.0 - L[i].sum()
    return L, B.stationary(q, Z=Z, beta=BETA), idx


ex, ey, ez = np.array([2.7, 0.0]), np.array([0.0, 2.7]), np.array([1.18, 1.18])
O = np.array([0.0, 0.0])
pos = lambda s: O + s[0] * ex + s[1] * ey + s[2] * ez
CENTROID = np.mean([pos(s) for s in B.STATES], axis=0)
LAB = {s: ('A' if s[0] == 0 else 'U') + ('D' if s[1] == 0 else 'C') + ('L' if s[2] == 0 else 'I')
       for s in B.STATES}
GOOD, BAD = (1, 1, 1), (1, 0, 0)
node_radius = lambda m: 0.18 + 0.55 * np.sqrt(m)


def draw(ax, q, title):
    L, pi, idx = chain(q)
    rad = {s: node_radius(pi[idx[s]]) for s in B.STATES}
    for s in B.STATES:
        for pop in range(3):
            t = list(s); t[pop] = 1 - t[pop]; t = tuple(t)
            if B.STATES.index(t) > idx[s]:
                p, r = pos(s), pos(t)
                ax.plot([p[0], r[0]], [p[1], r[1]], color='#cfd5db', lw=0.7, zorder=1)
    mx = max(L[i, j] for i in range(8) for j in range(8) if i != j)
    for s in B.STATES:
        i = idx[s]
        for pop in range(3):
            t = list(s); t[pop] = 1 - t[pop]; t = tuple(t)
            j = idx[t]
            if L[i, j] < 1.5e-4 or L[i, j] < L[j, i]:
                continue
            p, r = pos(s), pos(t)
            u = (r - p) / np.linalg.norm(r - p)
            a, b = p + u * (rad[s] + 0.03), r - u * (rad[t] + 0.05)
            frac = L[i, j] / mx
            strong = frac > 0.45
            ax.add_patch(FancyArrowPatch(a, b, arrowstyle='-|>', mutation_scale=6 + 5 * frac,
                                         lw=0.6 + 2.0 * frac, color=A_STR if strong else A_WEAK,
                                         alpha=0.95 if strong else 0.8, zorder=2,
                                         shrinkA=0, shrinkB=0))
    for s in B.STATES:
        p, m, R = pos(s), pi[idx[s]], rad[s]
        fc, ec, tc = C_GOOD if s == GOOD else C_BAD if s == BAD else C_MID
        ax.add_patch(Circle(p, R, facecolor=fc, edgecolor=ec, lw=0.9, zorder=3))
        if m >= 0.06:
            ax.text(p[0], p[1], f'{m*100:.0f}', ha='center', va='center',
                    fontsize=8.5, fontweight='bold', color=tc, zorder=5)
        if s == (1, 1, 0):
            d = np.array([-0.6, 0.8])
        elif s == (0, 0, 1):
            d = np.array([0.6, -0.8])
        else:
            d = p - CENTROID
        d = d / (np.linalg.norm(d) + 1e-9)
        lp = p + d * (R + 0.22)
        ha = 'center' if abs(d[0]) < 0.35 else ('left' if d[0] > 0 else 'right')
        va = 'center' if abs(d[1]) < 0.35 else ('bottom' if d[1] > 0 else 'top')
        txt = LAB[s] + (f' {m*100:.0f}%' if 0 < m < 0.06 and m >= 0.005 else '')
        ax.text(lp[0], lp[1], txt, ha=ha, va=va, fontsize=8, color='#2b2f33', zorder=5,
                fontweight=('bold' if s in (GOOD, BAD) else 'normal'))
    ax.text(*(O + ex * 0.5 + np.array([0.0, -0.62])), 'User', fontsize=8, color='#7d858d',
            style='italic', ha='center', va='center')
    ax.text(*(O + ey * 0.5 + np.array([-0.80, 0.0])), 'Concessionaire', fontsize=8, color='#7d858d',
            style='italic', ha='center', va='center', rotation=90)
    sp = pos((1, 0, 0)) + ez * 0.5 + np.array([0.55, -0.32])
    ax.text(sp[0], sp[1], 'State', fontsize=8, color='#7d858d', style='italic',
            ha='center', va='center', rotation=45)
    ax.text(-1.25, 5.45, title, fontsize=9, fontweight='bold', color='#1f2429', ha='left', va='top')
    ax.set_xlim(-1.55, 5.35)
    ax.set_ylim(-1.45, 5.55)
    ax.set_aspect('equal')
    ax.axis('off')


mechs = [(dict(P, phi=0.0), '(A) Baseline / punishment'),
         (dict(P, phi=0.0, bfo=1.2), '(B) Inspector reward'),
         (dict(P, phi=1.0, bR5=0.0), '(C) Exit valve open'),
         (dict(P, phi=1.0, bR5=1.9), '(D) Exit valve closed')]
fig = plt.figure(figsize=(TW, 0.98 * TW))
for k, (q, tt) in enumerate(mechs):
    ax = fig.add_axes([0.0 + 0.5 * (k % 2), 0.5 - 0.5 * (k // 2), 0.5, 0.5])
    draw(ax, q, tt)
save(fig, 'fig6_cube.png')
for q, tt in mechs:
    pi = B.stationary(q, Z=Z, beta=BETA)
    print(f'{tt:28s} pi(target)={pi[B.STATES.index(GOOD)]:.3f}  pi(UDL)={pi[B.STATES.index(BAD)]:.3f}')
