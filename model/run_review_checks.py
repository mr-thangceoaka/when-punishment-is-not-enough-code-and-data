"""Numerical robustness checks of the finite population and basin results.

  R1a  general fixation formula vs geometric formula (step rule)
  R1b  smooth logistic rule: switching point over y0 and steepness s
  R1c  reward size relative to the penalty (r = k*u)
  R1d  deterministic basins of the target state, punishment vs the state dependent rule
  R2a  the State keeps a share kappa of the fines: switching point and phi**
  R6   unequal population sizes Z_U, Z_C, Z_S
Writes review_checks.txt (console copy, read by figures/si_figs.py) and draft figures
figP_rule_robustness.png and figQ_basins.png.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import sys, io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import botppp_model as B
import fixation_general as G

out = io.StringIO()
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')

P = dict(B.DEFAULT, bR5=0.0)
phis = np.linspace(0, 1, 101)

log('=== R1a: general vs geometric fixation, step rule y0=1/2 ===')
for Z, beta in [(25, 0.1), (50, 0.5), (50, 10.0), (200, 2.0), (500, 10.0)]:
    old = np.array([B.pi_good(dict(P, scheme='adaptive', phi=ph), Z=Z, beta=beta) for ph in phis])
    new = np.array([G.pi_good(dict(P, phi=ph), Z=Z, beta=beta, rule='step') for ph in phis])
    log(f'  Z={Z:4d} beta={beta:5.1f}: max|diff| = {np.abs(old-new).max():.1e}')
log('  reason: with r_m=u_m, r_h=u_h the concessionaire payoff difference does not depend on w')

log('\n=== R1b: logistic rule w=1/(1+exp(-(y-y0)/s)), general formula, Z=50 beta=0.5 ===')
lo_p, _ = G.liftoff(P, phis, Z=50, beta=0.5, rule='punish')
lo_r, _ = G.liftoff(P, phis, Z=50, beta=0.5, rule='reward')
lo_a, _ = G.liftoff(P, phis, Z=50, beta=0.5, rule='step')
log(f'  punish {lo_p:.2f}  reward {lo_r:.2f}  step(y0=0.5) {lo_a:.2f}')
Y0 = [0.1, 0.25, 0.5, 0.75, 0.9]
S = [0.01, 0.05, 0.1, 0.2]
TAB_b = np.zeros((len(Y0), len(S)))
for i, y0 in enumerate(Y0):
    row = []
    for j, s in enumerate(S):
        TAB_b[i, j], _ = G.liftoff(P, phis, Z=50, beta=0.5, rule='logistic', y0=y0, s=s)
    log(f'  y0={y0:4.2f}: ' + '  '.join(f's={s}:{TAB_b[i,j]:.2f}' for j, s in enumerate(S)))
log('  same grid, beta=10 (strong selection):')
TAB_b10 = np.zeros_like(TAB_b)
for i, y0 in enumerate(Y0):
    for j, s in enumerate(S):
        TAB_b10[i, j], _ = G.liftoff(P, phis, Z=50, beta=10.0, rule='logistic', y0=y0, s=s)
    log(f'  y0={y0:4.2f}: ' + '  '.join(f's={s}:{TAB_b10[i,j]:.2f}' for j, s in enumerate(S)))

log('\n=== R1c: reward r_m+r_h = k*(u_m+u_h); condition r_m+r_h > c_P+g_m+g_h = 1.1 <=> k > 0.6875 ===')
K = [0.4, 0.6, 0.7, 0.8, 1.0, 1.25, 1.5]
TAB_c = {}
for k in K:
    q = dict(P, rm=0.8*k, rh=0.8*k)
    res = []
    for Z, beta in [(50, 0.5), (50, 10.0)]:
        l, _ = G.liftoff(q, phis, Z=Z, beta=beta, rule='step')
        res.append(l)
    TAB_c[k] = res
    log(f'  k={k:4.2f} (r_m+r_h={1.6*k:.2f}): lift-off Z=50 beta=0.5 -> {res[0]:.2f}; beta=10 -> {res[1]:.2f}')

log('\n=== R1d: deterministic basin of (1,1,1): share of 8^3 interior starts reaching it ===')
g = (np.arange(8) + 0.5)/8
X0 = np.array([[a, b, c] for a in g for b in g for c in g])
def basin(scheme, phi, rm=0.8, s=0.05, y0=0.5, T=1500.0):
    """Integrates the replicator system with the smooth weight; returns share ending near (1,1,1)."""
    S0 = X0.copy()
    q = dict(P, phi=phi, rm=rm, rh=rm)
    def f(Sx):
        x, y, z = Sx[:, 0], Sx[:, 1], Sx[:, 2]
        if scheme == 'punish': w = np.ones_like(y)
        else: w = 1.0/(1.0 + np.exp(-(y - y0)/s))
        m = 1 - phi*(1 - z)
        dU = q['bU']*m*(y + q['theta']*(1 - y))
        dC = (-(q['cP'] + q['gm']) - q['gh']*x*m
              + w*(q['um']*z + q['uh']*x*z) + (1 - w)*(q['rm']*z + q['rh']*x*z))
        dR = (-q['cR'] + phi*(q['bR'] - q['bR5'])*x
              - w*(q['v'] + q['lam'])*x*(1 - y)
              - (1 - w)*y*(q['rm'] + q['rh']*x))
        return np.stack([x*(1-x)*dU, y*(1-y)*dC, z*(1-z)*dR], axis=1)
    dt = 0.1
    for _ in range(int(T/dt)):
        k1 = f(S0); k2 = f(S0 + dt/2*k1); k3 = f(S0 + dt/2*k2); k4 = f(S0 + dt*k3)
        S0 = np.clip(S0 + dt/6*(k1 + 2*k2 + 2*k3 + k4), 0, 1)
    return np.mean(np.all(S0 > 0.95, axis=1))

rng = np.random.default_rng(1)
mx = 0.0
for _ in range(200):
    pt = rng.random(3); ph = rng.random()
    q = dict(P, phi=ph, scheme='adaptive')
    d_ref = B.deltas(q, pt)
    x, y, z = pt; w = B.instrument_weight(q, y)
    m = 1 - ph*(1 - z)
    dC = (-(q['cP']+q['gm']) - q['gh']*x*m + w*(q['um']*z + q['uh']*x*z) + (1-w)*(q['rm']*z + q['rh']*x*z))
    dR = -q['cR'] + ph*q['bR']*x - w*(q['v']+q['lam'])*x*(1-y) - (1-w)*y*(q['rm'] + q['rh']*x)
    mx = max(mx, abs(dC - d_ref[1]), abs(dR - d_ref[2]))
log(f'  closed-form adaptive gradients vs four-cell sums: max deviation {mx:.1e}')

PHB = np.round(np.arange(0.0, 1.0001, 0.05), 2)
BAS = {}
for sch in ('punish', 'adaptive'):
    BAS[sch] = np.array([basin(sch, ph) for ph in PHB])
    log(f'  {sch:8s}: ' + ' '.join(f'{ph:.2f}:{b:.2f}' for ph, b in zip(PHB, BAS[sch])))
BAS['adaptive_k06'] = np.array([basin('adaptive', ph, rm=0.8*0.6) for ph in PHB])
log('  adaptive, k=0.6 (reward too small): ' + ' '.join(f'{ph:.2f}:{b:.2f}' for ph, b in zip(PHB, BAS['adaptive_k06'])))

log('\n=== R2a: State keeps a share kappa of fines (row 3 + kappa(u_m+u_h), row 7 + kappa u_m) ===')
def T_kappa(p, w, kappa):
    T = B.payoff_cells(p, w=w)
    q = {**B.DEFAULT, **p}
    T = dict(T)
    u, c, s = T[(1, 0, 1)]; T[(1, 0, 1)] = (u, c, s + kappa*w*(q['um'] + q['uh']))
    u, c, s = T[(0, 0, 1)]; T[(0, 0, 1)] = (u, c, s + kappa*w*q['um'])
    return T
def pi_good_kappa(p, kappa, Z=50, beta=0.5):
    idx = {st: i for i, st in enumerate(B.STATES)}
    L = np.zeros((8, 8))
    T = T_kappa(p, 1.0, kappa)
    for st in B.STATES:
        i = idx[st]
        for pop in range(3):
            t = list(st); t[pop] = 1 - t[pop]; t = tuple(t)
            L[i, idx[t]] = B.fixation(T[t][pop] - T[st][pop], Z, beta)/3.0
        L[i, i] = 1 - L[i].sum()
    ev, V = np.linalg.eig(L.T); k = np.argmin(np.abs(ev - 1))
    pi = np.abs(np.real(V[:, k])); pi /= pi.sum()
    return pi[idx[(1, 1, 1)]]
for kappa in (0.0, 0.1, 0.25, 0.5, 1.0):
    phi2 = (P['cR'] + P['v'] + P['lam'] - kappa*(P['um'] + P['uh']))/P['bR']
    lo = []
    for beta in (0.5, 10.0):
        pg = np.array([pi_good_kappa(dict(P, phi=ph), kappa, beta=beta) for ph in phis])
        lo.append(phis[np.argmax(pg > 0.5)] if pg.max() > 0.5 else np.nan)
    log(f'  kappa={kappa:4.2f}: phi** = {phi2:+.3f}; lift-off beta=0.5 {lo[0]:.2f}, beta=10 {lo[1]:.2f}')
log('  gap closes when kappa(u_m+u_h) >= v+lambda, i.e. kappa >= 0.5 at canonical values')

log('\n=== R6: unequal population sizes (punish vs adaptive step), beta=0.5 and beta=10 ===')
def pi_good_Z(p, Zs, beta, rule):
    idx = {st: i for i, st in enumerate(B.STATES)}
    L = np.zeros((8, 8))
    for st in B.STATES:
        i = idx[st]; x, y, z = st
        for pop in range(3):
            t = list(st); t[pop] = 1 - t[pop]; t = tuple(t)
            Zp = Zs[pop]
            if pop == 1:
                mut = t[1]; dfs = []
                for j in range(1, Zp):
                    yy = j/Zp if mut == 1 else 1 - j/Zp
                    dfs.append(G.conc_delta(p, x, z, yy, mut, G.weight(rule, yy)))
                r = G.rho_general(dfs, beta)
            else:
                T = B.payoff_cells(p, w=G.weight(rule, float(y)))
                r = B.fixation(T[t][pop] - T[st][pop], Zp, beta)
            L[i, idx[t]] = r/3.0
        L[i, i] = 1 - L[i].sum()
    ev, V = np.linalg.eig(L.T); k = np.argmin(np.abs(ev - 1))
    pi = np.abs(np.real(V[:, k])); pi /= pi.sum()
    return pi[idx[(1, 1, 1)]]
ZSETS = [(50, 50, 50), (1000, 50, 50), (1000, 10, 20), (1000, 100, 10), (200, 5, 5), (10, 100, 100)]
for Zs in ZSETS:
    res = []
    for beta in (0.5, 10.0):
        for rule in ('punish', 'step'):
            pg = np.array([pi_good_Z(dict(P, phi=ph), Zs, beta, rule) for ph in phis])
            res.append(phis[np.argmax(pg > 0.5)] if pg.max() > 0.5 else np.nan)
    log(f'  (Z_U,Z_C,Z_S)={Zs}: beta=0.5 punish {res[0]:.2f} adaptive {res[1]:.2f} | '
        f'beta=10 punish {res[2]:.2f} adaptive {res[3]:.2f}')

plt.rcParams.update({'font.size': 12, 'axes.labelsize': 14, 'axes.titlesize': 13.5})
fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.6), gridspec_kw={'wspace': 0.32})
for ax, M, ttl in zip(axes, (TAB_b, TAB_b10), ('$Z=50$, $\\beta=0.5$', '$Z=50$, $\\beta=10$')):
    im = ax.imshow(M, cmap='viridis', vmin=0, vmax=0.6, aspect='auto')
    ax.set_xticks(range(len(S))); ax.set_xticklabels(S)
    ax.set_yticks(range(len(Y0))); ax.set_yticklabels(Y0)
    ax.set_xlabel('steepness $s$'); ax.set_ylabel('switching point $y_0$')
    ax.set_title(ttl)
    for i in range(len(Y0)):
        for j in range(len(S)):
            ax.text(j, i, f'{M[i,j]:.2f}', ha='center', va='center', color='w', fontsize=12, fontweight='bold')
cb = fig.colorbar(im, ax=axes, fraction=0.03); cb.set_label('lift off value of $\\phi$')
fig.savefig('figP_rule_robustness.png', dpi=200, bbox_inches='tight'); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.4, 4.4))
ax.plot(PHB, BAS['punish'], 'o-', color='tab:red', lw=2.2, label='punishment')
ax.plot(PHB, BAS['adaptive'], 's-', color='tab:green', lw=2.2, label='state dependent rule')
ax.plot(PHB, BAS['adaptive_k06'], '^:', color='tab:olive', lw=2.0,
        label='state dependent rule, reward $0.6\\times$ penalty')
ax.axvline(0.1, color='gray', ls=':', lw=1.6); ax.axvline(0.5, color='gray', ls='--', lw=1.6)
ax.set_xlabel('$\\phi$'); ax.set_ylabel('share of initial states reaching the target')
ax.set_ylim(-0.03, 1.03); ax.legend(fontsize=11, loc='center right')
fig.tight_layout(); fig.savefig('figQ_basins.png', dpi=200); plt.close(fig)

open('review_checks.txt', 'w', encoding='utf-8').write(out.getvalue())
log('\nwrote review_checks.txt, figP_rule_robustness.png, figQ_basins.png')
