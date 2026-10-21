"""Experiments for the two extensions in botppp_ext.py (SI, Sections 5 and 6).

E1  finite populations (rare mutations, extended precision): long run shares of the target (1,1,1)
    and of (Use, Comply, Lax) over (phi_o, phi_q), for punishment and the state dependent rule
E2  replicator dynamics: share of initial states ending with compliance (y > 0.95), and among them
    with audit, over (phi_o, phi_q), under punishment
E3  welfare of the eight pure states at the illustrative values; analytic differences
E4  discounted welfare along replicator paths from a grid of initial states, for no audit,
    punishment and the state dependent rule, over phi_o (phi_q = 0), with reward outlays
Writes extensions.txt, ext_E1.npy, ext_E2.npy, ext_E4.npy and draft figures figR_users.png and
figS_welfare.png.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import botppp_model as B
import botppp_ext as X
from run_grid_mp import rho_mp, gth, STATES, IDX

mp.mp.dps = 60
P = dict(B.DEFAULT, bR5=0.0)
out = []
def log(t): print(t, flush=True); out.append(t)

def weight(rule, y, y0=1/3, s=0.02):
    if rule == 'punish': return np.ones_like(np.asarray(y, float)) if np.ndim(y) else 1.0
    if rule == 'none': return np.ones_like(np.asarray(y, float)) if np.ndim(y) else 1.0
    if rule == 'step': return np.where(np.asarray(y) >= 0.5, 1.0, 0.0) if np.ndim(y) else (1.0 if y >= 0.5 else 0.0)
    return 1.0/(1.0 + np.exp(-(np.asarray(y) - y0)/s))

def pi_states(p, Z=50, beta=0.5, rule='punish'):
    L = mp.zeros(8, 8)
    for st in STATES:
        i = IDX[st]; x, y, z = st
        for pop in range(3):
            t = list(st); t[pop] = 1 - t[pop]; t = tuple(t)
            if pop == 1:
                mut = t[1]; dfs = []
                for j in range(1, Z):
                    yy = j/Z if mut == 1 else 1 - j/Z
                    T = X.payoff_cells_ext(p, w=float(weight(rule, yy)))
                    pc, pd = T[(x, 1, z)][1], T[(x, 0, z)][1]
                    dfs.append(pc - pd if mut == 1 else pd - pc)
            else:
                T = X.payoff_cells_ext(p, w=float(weight(rule, float(y))))
                dfs = [T[t][pop] - T[st][pop]]*(Z - 1)
            L[i, IDX[t]] = rho_mp(dfs, beta)/3
    return [float(v) for v in gth(L)]

PO = np.round(np.arange(0.0, 0.81, 0.05), 2)
PQ = np.round(np.arange(0.0, 0.81, 0.05), 2)
E1 = {}
for rule in ('punish', 'step'):
    A = np.full((len(PO), len(PQ), 2), np.nan)
    for i, po in enumerate(PO):
        for j, pq in enumerate(PQ):
            if po + pq > 1: continue
            pi = pi_states(dict(P, phi=float(po), phiq=float(pq)), rule=rule)
            A[i, j, 0] = pi[IDX[(1, 1, 1)]]; A[i, j, 1] = pi[IDX[(1, 1, 0)]]
    E1[rule] = A
np.save('ext_E1.npy', E1, allow_pickle=True)
log('=== E1: long run shares (Z=50, beta=0.5) ===')
for rule in ('punish', 'step'):
    A = E1[rule]
    log(f'-- {rule}: pi(1,1,1) | pi(1,1,0) at selected (phi_o, phi_q)')
    for po in (0.0, 0.05, 0.2, 0.4):
        i = list(PO).index(po)
        log(f'  phi_o={po:.2f}: ' + ' '.join(f'{PQ[j]:.2f}:{A[i,j,0]:.2f}|{A[i,j,1]:.2f}'
                                          for j in range(0, len(PQ), 2) if np.isfinite(A[i, j, 0])))
q_c = lambda po: (P['cP'] + P['gm'] + P['gh']*(1 - po))/(P['bP'] + P['gh'])
log(f'analytic phi_q^c(phi_o) = (c_P+g_m+g_h(1-phi_o))/(b_P+g_h): phi_o=0 {q_c(0):.3f}, 0.2 {q_c(0.2):.3f}, 0.4 {q_c(0.4):.3f}')

g = (np.arange(8) + 0.5)/8
X0 = np.array([[a, b, c] for a in g for b in g for c in g])

def field(S, q, rule):
    x, y, z = S[:, 0], S[:, 1], S[:, 2]
    w = weight(rule, y)
    po, pq = q['phi'], q['phiq']
    mo = 1 - po*(1 - z)
    dU = q['bU']*(y*mo + q['theta']*(1 - y)*(mo - pq))
    dC = (-(q['cP'] + q['gm']) + x*q['bP']*pq - x*q['gh']*(mo - pq)
          + w*(q['um']*z + q['uh']*x*z*(1 - pq)) + (1 - w)*(q['rm']*z + q['rh']*x*z))
    if rule == 'none':
        dR = -np.ones_like(z)*10.0
    else:
        dR = (-q['cR'] + po*(q['bR'] - q['bR5'])*x - w*(q['v'] + q['lam'])*x*(1 - y)
              - (1 - w)*y*(q['rm'] + q['rh']*x))
    return np.stack([x*(1 - x)*dU, y*(1 - y)*dC, z*(1 - z)*dR], axis=1), w

rng = np.random.default_rng(5); mx = 0.0
for _ in range(300):
    po = rng.random()*0.6; pq = rng.random()*(1 - po); S = rng.random((1, 3))
    q = {**X.EXT_DEFAULT, **dict(P, phi=po, phiq=pq)}
    for rule in ('punish', 'logistic'):
        F, w = field(S, q, rule)
        d = X.gradients(q, *S[0], w=float(w[0]))
        ref = np.array([S[0, 0]*(1 - S[0, 0])*d[0], S[0, 1]*(1 - S[0, 1])*d[1], S[0, 2]*(1 - S[0, 2])*d[2]])
        mx = max(mx, np.abs(F[0] - ref).max())
log(f'vectorised field vs four cell sums: max deviation {mx:.1e}')

def integrate(q, rule, S0, T=1500.0, dt=0.1, rho=None):
    """RK4; if rho is given also returns discounted welfare and reward outlay per start."""
    S = S0.copy(); n = int(T/dt)
    if rule == 'none':
        S[:, 2] = 0.0
    Wd = np.zeros(len(S)); Rd = np.zeros(len(S))
    for k in range(n):
        if rho is not None:
            _, w = field(S, q, rule)
            Wt, Rt = X.welfare(q, S[:, 0], S[:, 1], S[:, 2], w=w)
            disc = np.exp(-rho*k*dt)*dt
            Wd += disc*Wt; Rd += disc*Rt
        k1, _ = field(S, q, rule); k2, _ = field(S + dt/2*k1, q, rule)
        k3, _ = field(S + dt/2*k2, q, rule); k4, _ = field(S + dt*k3, q, rule)
        S = np.clip(S + dt/6*(k1 + 2*k2 + 2*k3 + k4), 0, 1)
    if rho is not None:
        return S, rho*Wd, rho*Rd
    return S

POb = np.round(np.arange(0.0, 0.61, 0.1), 2)
PQb = np.round(np.arange(0.0, 0.91, 0.1), 2)
E2 = np.full((len(POb), len(PQb), 2), np.nan)
for i, po in enumerate(POb):
    for j, pq in enumerate(PQb):
        if po + pq > 1: continue
        q = {**X.EXT_DEFAULT, **dict(P, phi=float(po), phiq=float(pq))}
        S = integrate(q, 'punish', X0)
        comp = S[:, 1] > 0.95
        E2[i, j, 0] = comp.mean(); E2[i, j, 1] = (comp & (S[:, 2] > 0.95)).mean()
np.save('ext_E2.npy', {'po': POb, 'pq': PQb, 'E2': E2}, allow_pickle=True)
log('=== E2: replicator, punishment: share of starts ending with compliance | with compliance and audit ===')
for i, po in enumerate(POb):
    log(f'  phi_o={po:.1f}: ' + ' '.join(f'{PQb[j]:.1f}:{E2[i,j,0]:.2f}|{E2[i,j,1]:.2f}'
                                        for j in range(len(PQb)) if np.isfinite(E2[i, j, 0])))

log('=== E3: welfare of the pure states (phi_o=0.2, phi_q=0; Ds=2, dq=0.7, kc=0.5, mu=0.3, eta=0.5) ===')
q = {**X.EXT_DEFAULT, **dict(P, phi=0.2)}
Wp = {}
for st in STATES:
    Wp[st] = float(X.welfare(q, *map(float, st), w=1.0)[0])
for st in sorted(Wp, key=lambda s: -Wp[s]):
    log(f'  {st}: W = {Wp[st]:+.3f}')

log('=== E4: discounted welfare (rho=0.01) averaged over 512 starts; phi_q=0 ===')
PH4 = np.round(np.arange(0.05, 0.651, 0.05), 2)
E4 = {r: np.zeros((len(PH4), 3)) for r in ('none', 'punish', 'logistic')}
for i, po in enumerate(PH4):
    q = {**X.EXT_DEFAULT, **dict(P, phi=float(po))}
    for rule in ('none', 'punish', 'logistic'):
        S, Wd, Rd = integrate(q, rule, X0, rho=0.01)
        E4[rule][i] = [Wd.mean(), Rd.mean(), np.mean(np.all(S > 0.95, axis=1))]
    log(f'  phi_o={po:.2f}: ' + '  '.join(f'{r}: W={E4[r][i,0]:+.3f} rew={E4[r][i,1]:.3f} target={E4[r][i,2]:.2f}'
                                        for r in ('none', 'punish', 'logistic')))
np.save('ext_E4.npy', {'phi': PH4, **E4}, allow_pickle=True)
log('=== E4 sensitivity at phi_o=0.25 ===')
for mu in (0.0, 0.3, 1.0):
    for rho in (0.005, 0.02):
        q = {**X.EXT_DEFAULT, **dict(P, phi=0.25, mu=mu)}
        r = {}
        for rule in ('punish', 'logistic'):
            _, Wd, Rd = integrate(q, rule, X0, rho=rho)
            r[rule] = (Wd.mean(), Rd.mean())
        log(f'  mu={mu:.1f} rho={rho:.3f}: punish W={r["punish"][0]:+.3f}; rule W={r["logistic"][0]:+.3f} (rewards {r["logistic"][1]:.3f})')

open('extensions.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')

plt.rcParams.update({'font.size': 12, 'axes.labelsize': 14, 'axes.titlesize': 13})
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), gridspec_kw={'wspace': 0.32})
for ax, (A, ttl) in zip(axes, [(E1['punish'][:, :, 0], 'punishment: target (Use, Comply, Audit)'),
                               (E1['punish'][:, :, 1], 'punishment: (Use, Comply, Lax)'),
                               (E1['step'][:, :, 0] + E1['step'][:, :, 1], 'state dependent rule: compliance, any audit')]):
    im = ax.imshow(A, origin='lower', cmap='viridis', vmin=0, vmax=1, aspect='auto',
                   extent=[PQ[0] - 0.025, PQ[-1] + 0.025, PO[0] - 0.025, PO[-1] + 0.025])
    pqv = np.linspace(0, 0.8, 100)
    ax.plot([q_c(po) for po in PO], PO, 'w--', lw=2)
    ax.axhline(0.1, color='w', ls=':', lw=1.8)
    ax.set_xlabel('$\\phi_q$ (users who observe performance)')
    ax.set_ylabel('$\\phi_o$ (users who observe oversight)')
    ax.set_title(ttl)
cb = fig.colorbar(im, ax=axes, fraction=0.02); cb.set_label('long run share')
fig.savefig('figR_users.png', dpi=200, bbox_inches='tight'); plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.4), gridspec_kw={'wspace': 0.28})
lab = {'none': 'no audit', 'punish': 'punishment', 'logistic': 'state dependent rule ($y_0=1/3$)'}
col = {'none': 'gray', 'punish': 'tab:red', 'logistic': 'tab:blue'}
for r in ('none', 'punish', 'logistic'):
    axes[0].plot(PH4, E4[r][:, 0], marker='o', color=col[r], lw=2.2, label=lab[r])
    axes[1].plot(PH4, E4[r][:, 1], marker='o', color=col[r], lw=2.2, label=lab[r])
for ax in axes:
    ax.axvline(0.1, color='gray', ls=':', lw=1.5); ax.axvline(0.5, color='gray', ls='--', lw=1.5)
    ax.set_xlabel('$\\phi_o$')
axes[0].set_ylabel('discounted welfare per period'); axes[0].legend(fontsize=10.5)
axes[1].set_ylabel('discounted reward outlay per period')
fig.savefig('figS_welfare.png', dpi=200, bbox_inches='tight'); plt.close(fig)
print('wrote extensions.txt figR_users.png figS_welfare.png')
