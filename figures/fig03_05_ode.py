"""SI figures fig3_reward.png, fig4_lending.png and fig5_regimes.png (replicator dynamics): the
inspector reward gives closed orbits, the effect of separating the lending role, and the three
regimes of the exit valve.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
from scipy.integrate import solve_ivp
from print_style import style, save, TW, plt
import botppp_model as B

style()
P = dict(B.DEFAULT)

p1 = dict(P, bfo=1.2, phi=0.0)
sol = solve_ivp(B.rhs, [0, 4000], [0.8, 0.8, 0.8], args=(p1,), rtol=1e-10, atol=1e-12,
                dense_output=True)
t = np.linspace(0, 4000, 40001)
Y = sol.sol(t)
fig = plt.figure(figsize=(TW, 2.05))
a0 = fig.add_axes([0.095, 0.20, 0.50, 0.66])
a1 = fig.add_axes([0.72, 0.20, 0.26, 0.66])
a0.plot(t, Y[1], lw=1.0, label='$y$ (Comply)')
a0.plot(t, Y[2], lw=1.0, label='$z$ (Audit)')
a0.set_xlabel('time')
a0.set_ylabel('frequency')
a0.set_xlim(0, 400)
a0.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False)
tail = t > 1500
a1.plot(Y[1][tail], Y[2][tail], lw=1.2)
a1.set_xlabel('$y$ (Comply)')
a1.set_ylabel('$z$ (Audit)')
save(fig, 'fig3_reward.png')
amp = Y[2][tail].max() - Y[2][tail].min()
print(f'[fig3_reward] amplitude of z (t>1500) = {amp:.3f}')

bfos = np.linspace(0, 2.0, 81)
audit_lvl = {}
for lam_, lab in [(P['lam'], 'State that also lends ($\\lambda=0.4$)'),
                  (0.0, 'lending role separated ($\\lambda=0$)')]:
    X0 = np.tile([0.8, 0.8, 0.8], (len(bfos), 1))
    q = {**P, 'lam': lam_, 'phi': 0.0, 'bfo': bfos}
    S = B.rk4_batch(X0, q, T=3000, dt=0.05)
    zbar = np.zeros(len(bfos))
    nwin = 4000
    Scur = S.copy()
    for k in range(nwin):
        Xc, Yc, Zc = Scur[:, 0], Scur[:, 1], Scur[:, 2]
        dU, dC, dR = B.deltas_vec(Xc, Yc, Zc, q)
        Scur = np.clip(Scur + 0.05 * np.stack([Xc * (1 - Xc) * dU, Yc * (1 - Yc) * dC,
                                               Zc * (1 - Zc) * dR], 1), 0, 1)
        zbar += Scur[:, 2]
    audit_lvl[lab] = zbar / nwin
fig = plt.figure(figsize=(0.8 * TW, 2.35))
ax = fig.add_axes([0.13, 0.18, 0.84, 0.60])
for lab, zb in audit_lvl.items():
    ax.plot(bfos, zb, lw=1.6, label=lab)
ax.axvline(P['v'] + P['lam'] + P['cR'], color='gray', lw=1.2, ls=':')
ax.axvline(P['v'] + P['cR'], color='gray', lw=1.2, ls='--')
ax.text(P['v'] + P['cR'] - 0.03, 0.55, '$v+c_R$', ha='right', va='center', color='0.3')
ax.text(P['v'] + P['lam'] + P['cR'] + 0.03, 0.55, '$v+\\lambda+c_R$', ha='left', va='center', color='0.3')
ax.set_xlabel('$b_{fo}$ (inspection reward)')
ax.set_ylabel('long run audit level $\\bar z$')
ax.set_xlim(0, 2)
ax.grid(alpha=0.25)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=1, frameon=False)
save(fig, 'fig4_lending.png')
lift = {k: bfos[np.argmax(v > 0.05)] for k, v in audit_lvl.items()}
print('[fig4_lending] lift off', lift)

runs = [("baseline ($\\phi=0$)", dict(P, phi=0.0), '#1f77b4', '-'),
        ("valve open ($\\phi=1$, $b_{R5}=0$)", dict(P, phi=1.0, bR5=0.0), '#ff7f0e', '-'),
        ("valve closed ($\\phi=1$, $b_{R5}=1.9>1.8$)", dict(P, phi=1.0, bR5=1.9), '#2ca02c', '--')]
fig = plt.figure(figsize=(0.8 * TW, 2.3))
ax = fig.add_axes([0.13, 0.18, 0.84, 0.58])
ends = []
for lab, q, col, ls in runs:
    s = solve_ivp(B.rhs, [0, 3000], [0.8, 0.8, 0.8], args=(q,), rtol=1e-10, atol=1e-12,
                  dense_output=True)
    tt = np.linspace(0, 40, 801)
    YY = s.sol(tt)
    ax.plot(tt, YY[1], label=lab, color=col, ls=ls, lw=1.8)
    ends.append((lab, np.round(s.sol(3000.0), 3)))
ax.set_xlabel('time')
ax.set_ylabel('$y$ (Comply)')
ax.set_xlim(0, 40)
ax.set_ylim(-0.03, 1.05)
ax.grid(alpha=0.25)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=1, frameon=False)
save(fig, 'fig5_regimes.png')
for lab, e in ends:
    print(f'[fig5_regimes] {lab}: end (x,y,z) = {e}')
