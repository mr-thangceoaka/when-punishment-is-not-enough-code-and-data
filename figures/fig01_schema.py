"""Figure 1 of the paper: (a) the three populations and how they are coupled; (b) how each
population revises its conduct. Drawn in inches at print size.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Polygon, Wedge, Circle
from print_style import style, save, TW, plt

style()

W, H = TW, 3.62
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis('off')
FT, FL = 8.5, 8


def icon_gov(cx, cy, c, k):
    R = lambda x, y, w, h: ax.add_patch(Rectangle((cx + k * x, cy + k * y), k * w, k * h, fc=c, ec='none', zorder=5))
    R(-0.80, -0.46, 1.60, 0.09); R(-0.68, -0.37, 1.36, 0.06)
    for i in range(5):
        R(-0.52 + i * 0.26 - 0.045, -0.31, 0.09, 0.33)
    R(-0.62, 0.02, 1.24, 0.075)
    ax.add_patch(Polygon([(cx - k * 0.70, cy + k * 0.095), (cx + k * 0.70, cy + k * 0.095), (cx, cy + k * 0.31)],
                         fc=c, ec='none', zorder=5))
    ax.add_patch(Wedge((cx, cy + k * 0.30), k * 0.185, 0, 180, fc=c, ec='none', zorder=5))


def icon_people(cx, cy, c, k):
    for dx, s in ((-0.40, 0.86), (0.40, 0.86), (0.0, 1.0)):
        ax.add_patch(Wedge((cx + k * dx, cy - k * 0.42 * s), k * 0.27 * s, 0, 180, fc=c, ec='white', lw=0.6, zorder=5))
        ax.add_patch(Circle((cx + k * dx, cy + k * 0.02 * s), k * 0.135 * s, fc=c, ec='white', lw=0.6, zorder=5))


def icon_towers(cx, cy, c, k):
    base = cy - k * 0.46
    for bx, bw, bh in [(-0.62, 0.34, 0.60), (-0.16, 0.32, 0.90), (0.28, 0.32, 0.48)]:
        ax.add_patch(Rectangle((cx + k * bx, base), k * bw, k * bh, fc=c, ec='none', zorder=5))
        for r in range(int(bh / 0.16)):
            for j in range(2):
                ax.add_patch(Rectangle((cx + k * (bx + 0.065 + j * 0.145), base + k * (0.07 + r * 0.16)),
                                       k * 0.055, k * 0.075, fc='white', ec='none', zorder=6))


def box(x, y, w, h, fc, ec, title, lines, icon=None, k=0.30):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0.03,rounding_size=0.08',
                                facecolor=fc, edgecolor=ec, linewidth=1.0, zorder=3))
    top = y + h / 2
    if icon is not None:
        icon(x, top - 0.20, ec, k)
        ty = top - 0.43
    else:
        ty = top - 0.11
    ax.text(x, ty, title, ha='center', va='center', fontsize=FT, fontweight='bold', zorder=4)
    for i, ln in enumerate(lines):
        ax.text(x, ty - 0.16 - 0.135 * i, ln, ha='center', va='center', fontsize=FL, zorder=4, color='0.15')


def arrow(p0, p1, color, rad=0.0, ls='-', lw=1.2):
    ax.add_patch(FancyArrowPatch(p0, p1, connectionstyle=f'arc3,rad={rad}', arrowstyle='-|>',
                                 mutation_scale=9, linewidth=lw, color=color, linestyle=ls, zorder=2))


def lab(x, y, text, color, ha='center', rot=0):
    ax.text(x, y, text, fontsize=FL, ha=ha, va='center', color=color, zorder=5, fontweight='bold',
            rotation=rot, linespacing=1.1,
            bbox=dict(boxstyle='round,pad=0.10', fc='white', ec='none', alpha=0.92))


C_S, C_U, C_C, C_R = '#b8860b', '#1d4ed8', '#15803d', '#6b7280'
SX, SY, SW, SH = 2.635, 2.93, 1.80, 0.88
UX, UY, UW, UH = 1.32, 1.42, 1.40, 1.06
CX, CY, CW, CH = 3.95, 1.42, 1.40, 1.06
RX, RY, RW, RH = 1.32, 0.25, 1.30, 0.40
MX = (UX + CX) / 2
box(SX, SY, SW, SH, '#fdf3d0', C_S, 'STATE', ['audits or stays lax', 'both regulator and lender'], icon_gov)
box(UX, UY, UW, UH, '#dbeafe', C_U, 'USERS', ['use the road or divert', 'some leave when', 'oversight is lax'],
    icon_people, k=0.28)
box(CX, CY, CW, CH, '#dcfce7', C_C, 'CONCESSIONAIRE', ['complies or cuts corners', 'gains from poor service'],
    icon_towers, k=0.32)
box(RX, RY, RW, RH, '#f3f4f6', C_R, 'Cheaper road', ['the way out for users'])

arrow((UX + UW / 2 + 0.04, UY + 0.18), (CX - CW / 2 - 0.04, CY + 0.18), C_U, rad=-0.25)
lab(MX, UY + 0.47, 'pay tolls', C_U)
arrow((CX - CW / 2 - 0.04, CY - 0.22), (UX + UW / 2 + 0.04, UY - 0.22), C_C, rad=-0.25)
lab(MX, UY - 0.50, 'service quality', C_C)
arrow((1.80, UY + UH / 2 + 0.04), (1.97, SY - SH / 2 - 0.03), C_U, rad=-0.1)
lab(2.05, 2.22, 'public value', C_U, ha='left')
arrow((3.30, SY - SH / 2 - 0.03), (3.47, CY + CH / 2 + 0.04), C_S, rad=-0.1)
lab(3.56, 2.24, 'audits and\npenalties', C_S, ha='left')
arrow((4.48, CY + CH / 2 + 0.04), (SX + SW / 2 + 0.03, 2.80), '#dc2626', rad=0.3, ls='--', lw=1.3)
lab(4.42, 3.27, 'sanctioning a\ndefector hurts\nthe State as its\nown lender', '#dc2626')
arrow((SX - SW / 2 - 0.03, 2.70), (1.05, UY + UH / 2 + 0.04), '#7c3aed', rad=0.25, ls=':', lw=1.3)
lab(0.86, 2.95, 'lax oversight\ndrives users\naway', '#7c3aed')
arrow((UX, UY - UH / 2 - 0.03), (RX, RY + RH / 2 + 0.03), C_R)
lab(UX + 0.10, (UY - UH / 2 + RY + RH / 2) / 2, 'divert', C_R, ha='left')
XL, YT, XE = 0.24, 3.53, SX - SW / 2 + 0.30
ax.plot([RX - RW / 2 - 0.03, XL, XL, XE], [RY, RY, YT, YT], color=C_R, lw=1.2, ls='--', zorder=2,
        solid_joinstyle='round', dash_joinstyle='round')
ax.add_patch(FancyArrowPatch((XE, YT), (XE, SY + SH / 2 + 0.035), arrowstyle='-|>', mutation_scale=9,
                             linewidth=1.2, color=C_R, zorder=2))
lab(0.10, 1.25, 'less value on the cheaper road', C_R, rot=90)
save(fig, 'fig1a_schema.png')

W2, H2 = TW, 1.95
fig = plt.figure(figsize=(W2, H2))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W2)
ax.set_ylim(0, H2)
ax.axis('off')
GREEN, RED, GREY = '#15803d', '#dc2626', '#6b7280'
YB = 0.50


def person(cx, cy, c, s=1.0):
    ax.add_patch(Wedge((cx, cy - 0.065 * s), 0.052 * s, 0, 180, fc=c, ec='white', lw=0.5, zorder=5))
    ax.add_patch(Circle((cx, cy + 0.006 * s), 0.028 * s, fc=c, ec='white', lw=0.5, zorder=5))


def bars(x0, y0, heights, color):
    for i, h in enumerate(heights):
        ax.add_patch(Rectangle((x0 + i * 0.15, y0), 0.115, h, fc=color, ec='white', lw=0.6, zorder=4))


ax.add_patch(FancyBboxPatch((1.09, 1.70), 2.80, 0.20, boxstyle='round,pad=0.03,rounding_size=0.06',
                            fc='#1f2937', ec='none', zorder=3))
ax.text(W2 / 2, 1.80, 'How each group changes its behaviour', ha='center', va='center', color='white',
        fontsize=9, fontweight='bold', zorder=4)
bars(0.42, YB, [0.18, 0.31, 0.47, 0.66], GREEN)
ax.add_patch(FancyArrowPatch((0.30, YB - 0.02), (0.30, YB + 0.78), arrowstyle='-|>', mutation_scale=10, lw=1.4,
                             color=GREEN, zorder=4))
ax.text(0.72, 1.45, 'Pays better than\nthe group average', ha='center', va='center', fontsize=8,
        fontweight='bold', color=GREEN, linespacing=1.15)
ax.text(0.72, 0.24, 'more and more\nadopt it', ha='center', va='center', fontsize=8, color='0.15',
        linespacing=1.15)
CYc = 0.95
ax.add_patch(Circle((W2 / 2, CYc), 0.47, fc='#f3f4f6', ec='0.6', lw=0.9, zorder=2))
rng = np.random.default_rng(7)
cols = [GREEN] * 9 + [GREY] * 7
pts = []
while len(pts) < len(cols):
    p = rng.uniform([-1.0, -1.0], [1.0, 1.0])
    if p[0] ** 2 + p[1] ** 2 <= 0.92:
        pts.append(p)
for (dx, dy), c in zip(pts, cols):
    person(W2 / 2 + 0.36 * dx, CYc + 0.36 * dy, c)
ax.text(W2 / 2, 0.24, 'a population that\ncopies what works', ha='center', va='center', fontsize=8,
        fontweight='bold', color='0.25', linespacing=1.15)
bars(W2 - 0.42 - 0.565, YB, [0.66, 0.47, 0.31, 0.18], RED)
ax.add_patch(FancyArrowPatch((W2 - 0.30, YB + 0.78), (W2 - 0.30, YB - 0.02), arrowstyle='-|>', mutation_scale=10,
                             lw=1.4, color=RED, zorder=4))
ax.text(W2 - 0.72, 1.45, 'Pays worse than\nthe group average', ha='center', va='center', fontsize=8,
        fontweight='bold', color=RED, linespacing=1.15)
ax.text(W2 - 0.72, 0.24, 'fewer and fewer\nkeep it', ha='center', va='center', fontsize=8, color='0.15',
        linespacing=1.15)
ax.add_patch(FancyArrowPatch((1.10, 1.05), (W2 / 2 - 0.52, 1.00), connectionstyle='arc3,rad=-0.25',
                             arrowstyle='-|>', mutation_scale=9, lw=1.1, color=GREY, zorder=3))
ax.add_patch(FancyArrowPatch((W2 / 2 + 0.52, 1.00), (W2 - 1.10, 1.05), connectionstyle='arc3,rad=-0.25',
                             arrowstyle='-|>', mutation_scale=9, lw=1.1, color=GREY, zorder=3))
save(fig, 'fig1b_replicator.png')
