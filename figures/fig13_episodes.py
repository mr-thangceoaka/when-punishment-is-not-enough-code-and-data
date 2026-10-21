"""SI figure fig13_episodes.png: the largest price concession offered in seven Vietnamese BOT
episodes of toll avoidance or protest (2015 to 2022), and how each episode ended.
Data: ../data/episodes.csv (sources in data/SOURCES.md and in the SI).

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import os
import numpy as np
from matplotlib.lines import Line2D
from print_style import style, save, TW, DATA, plt

style()

rows = []
with open(os.path.join(DATA, 'episodes.csv'), encoding='utf-8') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        rows.append(line.rstrip('\n'))
hdr = rows[0].split(',')
data = [dict(zip(hdr, r.split(','))) for r in rows[1:]]
cuts = np.array([float(d['partial_cut_pct']) for d in data])
N = len(data)
order = sorted(range(N), key=lambda i: cuts[i])
res_color = {'full_exempt_local': 'C0', 'station_removal_buyout': 'C3', 'relocation': 'C4',
             'shutdown': 'C1', 'admin_ban': 'C2', 'enforcement_abandoned': 'C5'}
res_label = {'full_exempt_local': 'full local exemption',
             'station_removal_buyout': 'collection halted; State purchase proposed',
             'relocation': 'station relocated', 'shutdown': 'collection suspended',
             'admin_ban': 'administrative ban', 'enforcement_abandoned': 'enforcement abandoned'}
disp = {'CauRac': 'Cau Rac', 'HacTri': 'Hac Tri', 'BenThuy1': 'Ben Thuy 1', 'CaiLay': 'Cai Lay',
        'DT391_QL5': 'DT391/QL5*', 'TanDe': 'Tan De', 'T2_QL91': 'T2/QL91'}
fig = plt.figure(figsize=(TW, 2.55))
ax = fig.add_axes([0.16, 0.17, 0.82, 0.64])
for yy, i in enumerate(order):
    d = data[i]
    c = res_color[d['resolution']]
    ax.plot([0, cuts[i]], [yy, yy], '-', color='0.75', lw=0.9, zorder=1)
    ax.scatter([cuts[i]], [yy], s=40, color=c, zorder=3, marker='X' if cuts[i] > 0 else 'o',
               edgecolor='k', linewidth=0.5)
    ax.text(cuts[i] + 3, yy, res_label[d['resolution']], ha='left', va='center', fontsize=8, color=c)
ax.set_yticks(range(N))
ax.set_yticklabels([disp.get(data[i]['episode'], data[i]['episode']) for i in order])
ax.tick_params(axis='y', length=0)
ax.set_xlim(-3, 150)
ax.set_ylim(-0.6, N - 0.4)
ax.set_xlabel("largest price concession offered (%)")
ax.grid(axis='x', alpha=0.3)
leg = [Line2D([0], [0], marker='X', color='w', markerfacecolor='0.6', markeredgecolor='k', markersize=7,
              label='concession offered; avoidance or protest continued'),
       Line2D([0], [0], marker='o', color='w', markerfacecolor='0.6', markeredgecolor='k', markersize=6,
              label='no concession offered')]
fig.legend(handles=leg, loc="upper left", bbox_to_anchor=(0.15, 1.0), ncol=1, frameon=False,
           columnspacing=1.2)
save(fig, 'fig13_episodes.png')
