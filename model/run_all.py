"""Base runs of the model at the illustrative values.

Computes and caches the data of two print figures: the share of initial states that reach the
target over (phi, b_R5) (figH_basin.npy, Figure 2 of the paper) and the long run share of the
target in finite populations over (phi, b_R5) (figK_pigood.npy). Prints checks of the analytic
thresholds. Also writes draft figures (fig1_reward.png to fig6_calibration.png) that are not used
in the paper.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import botppp_model as B

P = dict(B.DEFAULT)
GOOD = (1,1,1)
rng = np.random.default_rng(1)
print("Canonical params:", P)
print(f"Conditions: um+uh={P['um']+P['uh']} > cP+gm+gh={P['cP']+P['gm']+P['gh']}; "
      f"valve ceiling bR5*={P['bR']-P['cR']}; phi*(bR5=0)={P['cR']/P['bR']}; "
      f"phi**(bR5=0)={(P['cR']+P['v']+P['lam'])/P['bR']}")

p1 = dict(P, bfo=1.2, phi=0.0)
sol = solve_ivp(B.rhs, [0,4000], [0.8,0.8,0.8], args=(p1,), rtol=1e-10, atol=1e-12,
                dense_output=True)
t = np.linspace(0,4000,4000); Y = sol.sol(t)
fig, axes = plt.subplots(1, 2, figsize=(12.2,4.8))
axes[0].plot(t, Y[1], lw=1.4, label='y (Comply)'); axes[0].plot(t, Y[2], lw=1.4, label='z (Audit)')
axes[0].set_xlabel('time', fontsize=15); axes[0].set_ylabel('frequency', fontsize=15)
axes[0].tick_params(axis='both', labelsize=13); axes[0].legend(fontsize=13.5, framealpha=0.95)
axes[0].set_title(f"Reward variant: $b_{{fo}}$=1.2 > $v+\\lambda+c_R$={P['v']+P['lam']+P['cR']:.1f}", fontsize=15)
tail = t > 1500
axes[1].plot(Y[1][tail], Y[2][tail], lw=1.6)
axes[1].set_xlabel('y (Comply)', fontsize=15); axes[1].set_ylabel('z (Audit)', fontsize=15)
axes[1].tick_params(axis='both', labelsize=13); axes[1].set_title('closed orbit in the (y,z) plane', fontsize=15)
fig.tight_layout(); fig.savefig('fig1_reward.png', dpi=200); plt.close(fig)
amp = Y[2][tail].max()-Y[2][tail].min()
print(f"[Fig1] amplitude z (t>1500) = {amp:.3f} -> {'LIMIT CYCLE ok' if amp>0.05 else 'NO oscillation'}")

bfos = np.linspace(0, 2.0, 81)
audit_lvl = {}
for lam_, lab in [(P['lam'], 'creditor State ($\\lambda$=0.4)'), (0.0, 'lending separated ($\\lambda$=0)')]:
    X0 = np.tile([0.8,0.8,0.8], (len(bfos),1))
    q = {**P, 'lam':lam_, 'phi':0.0, 'bfo':bfos}
    S = B.rk4_batch(X0, q, T=3000, dt=0.05)
    zbar = np.zeros(len(bfos)); nwin = 4000
    Scur = S.copy()
    for k in range(nwin):
        Xc, Yc, Zc = Scur[:,0], Scur[:,1], Scur[:,2]
        dU,dC,dR = B.deltas_vec(Xc,Yc,Zc,q)
        Scur = np.clip(Scur + 0.05*np.stack([Xc*(1-Xc)*dU, Yc*(1-Yc)*dC, Zc*(1-Zc)*dR],1), 0, 1)
        zbar += Scur[:,2]
    audit_lvl[lab] = zbar/nwin
fig, ax = plt.subplots(figsize=(7.6,5.0))
for lab, zb in audit_lvl.items(): ax.plot(bfos, zb, lw=2.6, label=lab)
ax.axvline(P['v']+P['lam']+P['cR'], color='gray', lw=2, ls=':', label='$v+\\lambda+c_R$')
ax.axvline(P['v']+P['cR'], color='gray', lw=2, ls='--', label='$v+c_R$')
ax.set_xlabel('$b_{fo}$ (inspection reward)', fontsize=16); ax.set_ylabel('long-run audit level $\\bar z$', fontsize=16)
ax.tick_params(axis='both', labelsize=13); ax.grid(alpha=0.25)
ax.legend(fontsize=13, framealpha=0.95); ax.set_title('Creditor loss shifts the activation threshold by exactly $\\lambda$', fontsize=15)
fig.tight_layout(); fig.savefig('fig2_lending.png', dpi=200); plt.close(fig)
def liftoff(zb):
    j = np.argmax(zb > 0.05); return bfos[j] if zb.max() > 0.05 else np.nan
l1 = liftoff(audit_lvl[list(audit_lvl)[0]]); l2 = liftoff(audit_lvl[list(audit_lvl)[1]])
print(f"[Fig2] lift-off: creditor = {l1:.3f}, separated = {l2:.3f}, gap = {l1-l2:.3f} (expected lambda = {P['lam']})")

fig, ax = plt.subplots(figsize=(9.2,5.2))
runs = [("baseline ($\\phi$=0)", dict(P, phi=0.0), '#1f77b4', '-'),
        ("valve OPEN ($\\phi$=1, $b_{R5}$=0)", dict(P, phi=1.0, bR5=0.0), '#ff7f0e', '-'),
        ("valve CLOSED ($\\phi$=1, $b_{R5}$=1.9 > 1.8)", dict(P, phi=1.0, bR5=1.9), '#2ca02c', '--')]
ends = []
for lab, q, col, ls in runs:
    s = solve_ivp(B.rhs, [0,3000], [0.8,0.8,0.8], args=(q,), rtol=1e-10, atol=1e-12, dense_output=True)
    tt = np.linspace(0,3000,1500); YY = s.sol(tt)
    ax.plot(tt, YY[1], label=f'y, {lab}', color=col, ls=ls, lw=3.0)
    ends.append((lab, np.round(YY[:,-1],3)))
ax.set_xlabel('time', fontsize=17); ax.set_ylabel('y (Comply)', fontsize=17)
ax.tick_params(axis='both', labelsize=14); ax.set_ylim(-0.03, 1.05); ax.grid(alpha=0.25)
ax.legend(fontsize=13.5, loc='center right', framealpha=0.95)
ax.set_title('Three regimes from the same start (0.8, 0.8, 0.8)', fontsize=18)
fig.tight_layout(); fig.savefig('fig3_regimes.png', dpi=200); plt.close(fig)
for lab, e in ends: print(f"[Fig3] {lab}: end (x,y,z) = {e}")

from matplotlib.lines import Line2D
def _field(Y, Z, q):
    _, dC, dR = B.deltas_vec(1.0, Y, Z, q)
    return Y*(1-Y)*dC, Z*(1-Z)*dR
def _kind(y, z, q, h=1e-5):
    def F(a, b):
        _, dC, dR = B.deltas_vec(1.0, a, b, q); return np.array([a*(1-a)*dC, b*(1-b)*dR])
    J = np.column_stack([(F(y+h,z)-F(y-h,z))/(2*h), (F(y,z+h)-F(y,z-h))/(2*h)])
    ev = np.linalg.eigvals(J).real
    return 'stable' if np.all(ev < -1e-9) else 'unstable' if np.all(ev > 1e-9) else 'saddle'
_g = np.linspace(0, 1, 400); _Yg, _Zg = np.meshgrid(_g, _g)
fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.5), constrained_layout=True)
strm = None
for ax, (lab, q) in zip(axes, [(r"Valve open  ($b_{R5}=0$)", dict(P, phi=1.0, bR5=0.0)),
                               (r"Valve closed  ($b_{R5}=1.9$)", dict(P, phi=1.0, bR5=1.9))]):
    U, V = _field(_Yg, _Zg, q); spd = np.hypot(U, V)
    strm = ax.streamplot(_g, _g, U, V, color=spd, cmap='viridis',
                         density=0.95, linewidth=1.0, arrowsize=0.9,
                         broken_streamlines=False, zorder=1)
    zc = (q['cP']+q['gm'])/(q['um']+q['uh']-q['gh'])
    if 0 < zc < 1: ax.axhline(zc, ls=(0,(6,4)), color='crimson', lw=1.5, alpha=.95, zorder=5)
    fpkw = dict(zorder=7, clip_on=False)
    for fy, fz in [(0,0),(1,0),(0,1),(1,1)]:
        k = _kind(fy, fz, q)
        ax.scatter([fy],[fz], s=210, marker=('o' if k!='saddle' else 's'),
                   facecolor='white', edgecolor='none', **fpkw)
        if k == 'stable':
            halo = 'seagreen' if (fy, fz) == (1,1) else 'firebrick'
            ax.scatter([fy],[fz], s=560, facecolor=halo, alpha=.28, edgecolor='none',
                       zorder=6, clip_on=False)
            ax.scatter([fy],[fz], s=165, marker='o', facecolor='black',
                       edgecolor='white', linewidths=1.8, **fpkw)
        elif k == 'unstable':
            ax.scatter([fy],[fz], s=150, marker='o', facecolor='white',
                       edgecolor='black', linewidths=1.8, **fpkw)
        else:
            ax.scatter([fy],[fz], s=130, marker='s', facecolor='0.82',
                       edgecolor='black', linewidths=1.6, **fpkw)
    ax.annotate('honest\n$(1,1)$', (1,1), (0.80,0.865), fontsize=13, color='darkgreen', ha='center', zorder=7)
    ax.annotate('misconduct\n$(0,0)$', (0,0), (0.17,0.105), fontsize=13, color='firebrick', ha='center', zorder=7)
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_aspect('equal')
    ax.set_xlabel(r'$y$  (fraction Comply)', fontsize=15); ax.set_ylabel(r'$z$  (fraction Audit)', fontsize=15)
    ax.tick_params(axis='both', labelsize=12.5)
    ax.set_title(lab, fontsize=15, pad=8)
cb = fig.colorbar(strm.lines, ax=axes, shrink=0.86, pad=0.015)
cb.set_label(r'flow speed  $\|(\dot y,\dot z)\|$', fontsize=14); cb.ax.tick_params(labelsize=12)
fig.legend(handles=[
    Line2D([],[],marker='o',color='w',mfc='black',mec='white',ms=11,label='stable (attractor)'),
    Line2D([],[],marker='o',color='w',mfc='white',mec='black',ms=10,label='unstable'),
    Line2D([],[],marker='s',color='w',mfc='0.82',mec='black',ms=10,label='saddle'),
    Line2D([],[],ls=(0,(6,4)),color='crimson',label=r'$\dot y=0$ nullcline')],
    loc='outside lower center', ncol=4, fontsize=12.5, framealpha=.92)
fig.savefig('fig4_plane.png', dpi=200); plt.close(fig)
print("[Fig4] done")

nph, nbr = 26, 26
phis = np.linspace(0, 1, nph); brs = np.linspace(0, 2.0, nbr)
g = np.linspace(0.15, 0.85, 3)
starts27 = np.array([(a,b,c) for a in g for b in g for c in g])
PH, BR = np.meshgrid(phis, brs, indexing='ij')
Ncell = nph*nbr
X0 = np.tile(starts27, (Ncell,1))
qv = {**P,
      'phi': np.repeat(PH.ravel(), len(starts27)),
      'bR5': np.repeat(BR.ravel(), len(starts27))}
S = B.rk4_batch(X0, qv, T=2500, dt=0.1)
good = ((np.abs(S - np.array(GOOD)) < 0.05).all(axis=1)).reshape(Ncell, len(starts27))
basin = good.mean(axis=1).reshape(nph, nbr)
np.save('figH_basin.npy', basin)
fig, ax = plt.subplots(figsize=(7.2,5.2))
im = ax.imshow(basin.T, origin='lower', extent=[0,1,0,2], aspect='auto', cmap='viridis', vmin=0, vmax=1)
ph = np.linspace(P['cR']/P['bR']+1e-6, 1, 200)
ax.plot(ph, P['bR'] - P['cR']/ph, c='cyan', lw=3, label="analytic $b_{R5} = b_R - c_R/\\phi$  ($\\phi^*$)")
ax.set_xlabel('$\\phi$ (discerning fraction)', fontsize=17); ax.set_ylabel('$b_{R5}$ (escape-route toll)', fontsize=17)
ax.tick_params(axis='both', labelsize=13)
ax.set_ylim(0,2); ax.legend(loc='lower right', fontsize=13.5, framealpha=0.97)
cb = fig.colorbar(im); cb.set_label('basin volume of good state', fontsize=15); cb.ax.tick_params(labelsize=12)
ax.set_title('Infinite-population basin volume', fontsize=18, pad=14)
fig.tight_layout(); fig.savefig('figH_phi_heatmap.png', dpi=200); plt.close(fig)
err = []
for j, br in enumerate(brs):
    col = basin[:, j]
    if col.max() > 0:
        phi_num = phis[np.argmax(col > 0)]
        phi_ana = P['cR']/(P['bR']-br) if P['bR']-br > P['cR'] else np.inf
        if np.isfinite(phi_ana): err.append(abs(phi_num - phi_ana))
print(f"[FigH] boundary deviation, numeric vs analytic (mean over bR5) = {np.mean(err):.3f} (phi grid step {phis[1]-phis[0]:.3f})")

fig, ax = plt.subplots(figsize=(7.2,4.7))
nstart = 200
starts_r = rng.uniform(0.05, 0.95, size=(nstart,3))
phis_f = np.linspace(0, 1, 51)
for br, col in [(0.0,'tab:blue'), (1.1,'tab:orange'), (1.9,'tab:red')]:
    frac = np.zeros(len(phis_f))
    X0 = np.tile(starts_r, (len(phis_f),1))
    qv = {**P, 'bR5': br, 'phi': np.repeat(phis_f, nstart)}
    S = B.rk4_batch(X0, qv, T=2500, dt=0.1)
    good = ((np.abs(S - np.array(GOOD)) < 0.05).all(axis=1)).reshape(len(phis_f), nstart)
    frac = good.mean(axis=1)
    pstar = P['cR']/(P['bR']-br) if P['bR']-br > 0 else np.inf
    ax.plot(phis_f, frac, color=col, lw=2.8, label=f"$b_{{R5}}$={br}  ($\\phi^*$={pstar:.2f})")
    if pstar <= 1: ax.axvline(pstar, color=col, ls=':', lw=2, alpha=0.7)
ax.set_xlabel('$\\phi$', fontsize=17); ax.set_ylabel('P(good state | random start)', fontsize=16)
ax.tick_params(axis='both', labelsize=13); ax.grid(alpha=0.25)
ax.legend(fontsize=13.5, framealpha=0.95); ax.set_title('Rescue switches on at $\\phi^*$', fontsize=18)
fig.tight_layout(); fig.savefig('figI_phi_bifurcation.png', dpi=200); plt.close(fig)
print("[FigI] done")

mechs = [("baseline", dict(P, phi=0.0)),
         ("reward $b_{fo}$=1.2", dict(P, phi=0.0, bfo=1.2)),
         ("valve OPEN", dict(P, phi=1.0, bR5=0.0)),
         ("valve CLOSED", dict(P, phi=1.0, bR5=1.9))]
labels8 = ['ADL','ADA','ACL','ACA','UDL','UDA','UCL','UCA']
order = [(0,0,0),(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1)]
fig, axes = plt.subplots(1, 4, figsize=(15,3.6), sharey=True)
pis = {}
for ax,(lab,q) in zip(axes, mechs):
    pi = B.stationary(q, Z=50, beta=0.5)
    pis[lab] = pi
    vals = [pi[B.STATES.index(s)] for s in order]
    cols = ['tab:green' if s==GOOD else ('tab:red' if s[1]==0 and s[2]==0 else 'tab:blue') for s in order]
    ax.bar(range(8), vals, color=cols)
    ax.set_xticks(range(8)); ax.set_xticklabels(labels8, fontsize=7, rotation=45)
    ax.set_title(lab, fontsize=10)
axes[0].set_ylabel('stationary mass $\\pi$')
fig.suptitle('Finite-population stationary distribution (Z=50, $\\beta$=0.5); green = good state (Use, Comply, Audit)')
fig.tight_layout(); fig.savefig('fig5_stationary.png', dpi=130); plt.close(fig)
for lab, pi in pis.items():
    pg = pi[B.STATES.index(GOOD)]
    dl = sum(pi[B.STATES.index(s)] for s in B.STATES if s[1]==0 and s[2]==0)
    print(f"[Fig5] {lab:16s}: pi(good) = {pg:.3f}, pi(Defect&Lax) = {dl:.3f}")

nph2, nbr2 = 26, 26
phis2 = np.linspace(0,1,nph2); brs2 = np.linspace(0,2,nbr2)
PG = np.zeros((nph2,nbr2))
for i,ph_ in enumerate(phis2):
    for j,br_ in enumerate(brs2):
        PG[i,j] = B.pi_good(dict(P, phi=ph_, bR5=br_), Z=50, beta=0.5)
np.save('figK_pigood.npy', PG)
fig, ax = plt.subplots(figsize=(7.2,5.2))
im = ax.imshow(PG.T, origin='lower', extent=[0,1,0,2], aspect='auto', cmap='viridis', vmin=0, vmax=1)
ph = np.linspace(P['cR']/P['bR']+1e-6, 1, 300)
ax.plot(ph, P['bR'] - P['cR']/ph, c='cyan', lw=3, label='local stability  $\\phi^*$')
ph2 = np.linspace((P['cR']+P['v']+P['lam'])/P['bR']+1e-6, 1, 300)
ax.plot(ph2, P['bR'] - (P['cR']+P['v']+P['lam'])/ph2, color='k', ls='--', lw=3, label='reachability  $\\phi^{**}$')
ax.set_xlabel('$\\phi$', fontsize=17); ax.set_ylabel('$b_{R5}$', fontsize=17)
ax.tick_params(axis='both', labelsize=13); ax.set_ylim(0,2)
ax.legend(loc='lower right', fontsize=13.5, framealpha=0.97)
cb = fig.colorbar(im); cb.set_label('$\\pi$(good)', fontsize=15); cb.ax.tick_params(labelsize=12)
ax.set_title('Finite-population mass on the good state (Z=50, $\\beta$=0.5)', fontsize=15.5, pad=14)
fig.tight_layout(); fig.savefig('figK_phi_finite_heatmap.png', dpi=200); plt.close(fig)
print("[FigK] done")

fig, ax = plt.subplots(figsize=(7.2,4.7))
phis3 = np.linspace(0,1,101)
for br, col in [(0.0,'tab:blue'), (1.1,'tab:orange')]:
    pg = [B.pi_good(dict(P, phi=ph_, bR5=br), Z=50, beta=0.5) for ph_ in phis3]
    ps  = P['cR']/(P['bR']-br)
    pss = (P['cR']+P['v']+P['lam'])/(P['bR']-br)
    ax.plot(phis3, pg, color=col, lw=2.8, label=f"$b_{{R5}}$={br}")
    for val, ls in [(ps,':'),(pss,'--')]:
        if val <= 1: ax.axvline(val, color=col, ls=ls, lw=2, alpha=0.7)
ax.set_xlabel('$\\phi$', fontsize=17); ax.set_ylabel('$\\pi$(good)', fontsize=16)
ax.tick_params(axis='both', labelsize=13); ax.grid(alpha=0.25)
ax.legend(fontsize=13, loc='upper left', framealpha=0.95)
ax.set_title('$\\pi$(good) vs $\\phi$: lift-off strictly above $\\phi^*$', fontsize=16)
fig.tight_layout(); fig.savefig('figL_pi_vs_phi.png', dpi=200); plt.close(fig)
pg_arr = np.array([B.pi_good(dict(P, phi=ph_, bR5=0.0)) for ph_ in phis3])
lift = phis3[np.argmax(pg_arr > 0.1)]
print(f"[FigL] bR5=0: lift-off numeric = {lift:.2f}; phi* = 0.10, phi** = 0.50")

fig, ax = plt.subplots(figsize=(7,4.5))
r = np.linspace(0, 0.95, 200)
phi_s  = (P['cR']/P['bR'])/(1-r)
phi_ss = ((P['cR']+P['v']+P['lam'])/P['bR'])/(1-r)
ax.plot(r, np.minimum(phi_s,1.05), label='$\\phi^*$ (stability)')
ax.plot(r, np.minimum(phi_ss,1.05), label='$\\phi^{**}$ (reachability)')
ax.fill_between(r, np.minimum(phi_s,1.05), np.minimum(phi_ss,1.05), alpha=0.2, label='trapped band')
ax.axvline(0.55, color='k', ls='--', label='VIDIFI 2025: $b_{R5}/b_R \\approx 0.55$')
ax.axhline(1.0, color='gray', lw=0.5)
ax.set_xlabel('$b_{R5}/b_R$'); ax.set_ylabel('$\\phi$ threshold'); ax.set_ylim(0,1.08)
ax.legend(fontsize=8); ax.set_title('Thresholds vs toll ratio (illustrative $c_R, v+\\lambda$ in revenue units)')
fig.tight_layout(); fig.savefig('fig6_calibration.png', dpi=130); plt.close(fig)
r55s, r55ss = (P['cR']/P['bR'])/0.45, ((P['cR']+P['v']+P['lam'])/P['bR'])/0.45
print(f"[Fig6] at r=0.55: phi* = {r55s:.3f}, phi** = {r55ss:.3f} ({'TRAPPED (phi**>1)' if r55ss>1 else 'assembleable'})")
print("ALL DONE")
