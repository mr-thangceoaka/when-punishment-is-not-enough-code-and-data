"""Deterministic basins of the target state under the state dependent rule.

Share of initial states that reach the target under the replicator dynamics, as a function of phi
and of the switching point y0. Along the face x = 1, z = 1 the State's fitness difference under a
rule w(y) is Delta_R(y) = phi (b_R - b_R5) - c_R - c(y), with
c(y) = w(y) (v + lambda)(1 - y) + (1 - w(y)) y (r_m + r_h): under punishment the State pays
v + lambda per defector found, under reward r_m + r_h per complier found. Audit is favoured along
the whole path if and only if phi (b_R - b_R5) > c_R + max_y c(y). For a step rule the maximum is
smallest at y0 = (v + lambda)/(v + lambda + r_m + r_h), which gives the threshold phi_dagger.
Writes basins.npy, basins.txt and a draft figure figQ_basins.png.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import botppp_model as B

P = dict(B.DEFAULT, bR5=0.0)
L = P['v'] + P['lam']; R = P['rm'] + P['rh']
y0_opt = L/(L + R)

def phi_path(y0, rsum=R):
    return (P['cR'] + max(y0*rsum, (1 - y0)*L))/(P['bR'] - P['bR5'])

g = (np.arange(10) + 0.5)/10
X0 = np.array([[a, b, c] for a in g for b in g for c in g])

def basin(phi, rule, y0=0.5, s=0.02, T=3000.0, dt=0.1, rm=None):
    q = dict(P, phi=phi)
    if rm is not None: q.update(rm=rm, rh=rm)
    S = X0.copy()
    def f(Sx):
        x, y, z = Sx[:, 0], Sx[:, 1], Sx[:, 2]
        w = np.ones_like(y) if rule == 'punish' else 1.0/(1.0 + np.exp(-(y - y0)/s))
        m = 1 - phi*(1 - z)
        dU = q['bU']*m*(y + q['theta']*(1 - y))
        dC = (-(q['cP'] + q['gm']) - q['gh']*x*m + w*(q['um']*z + q['uh']*x*z)
              + (1 - w)*(q['rm']*z + q['rh']*x*z))
        dR = (-q['cR'] + phi*(q['bR'] - q['bR5'])*x - w*(q['v'] + q['lam'])*x*(1 - y)
              - (1 - w)*y*(q['rm'] + q['rh']*x))
        return np.stack([x*(1 - x)*dU, y*(1 - y)*dC, z*(1 - z)*dR], axis=1)
    for _ in range(int(T/dt)):
        k1 = f(S); k2 = f(S + dt/2*k1); k3 = f(S + dt/2*k2); k4 = f(S + dt*k3)
        S = np.clip(S + dt/6*(k1 + 2*k2 + 2*k3 + k4), 0, 1)
    return float(np.mean(np.all(S > 0.95, axis=1)))

if __name__ == '__main__':
    lines = []
    def log(t): print(t); lines.append(t)
    log(f'phi* = {P["cR"]/P["bR"]:.3f}, phi** = {(P["cR"]+L)/P["bR"]:.3f}')
    log(f'path condition: y0=0.5 -> phi_path = {phi_path(0.5):.3f}; '
        f'best y0 = {y0_opt:.3f} -> phi_path = {phi_path(y0_opt):.3f}')
    PH = np.round(np.arange(0.05, 0.751, 0.05), 2)
    cases = [('punish', None), ('step', 0.5), ('step', round(y0_opt, 3)), ('step', 0.2)]
    res = {}
    for rule, y0 in cases:
        key = rule if y0 is None else f'y0={y0}'
        res[key] = np.array([basin(ph, rule, y0 if y0 else 0.5) for ph in PH])
        log(f'{key:10s}: ' + ' '.join(f'{p:.2f}:{b:.2f}' for p, b in zip(PH, res[key])))
        if y0 is not None:
            log(f'{"":10s}  path condition phi > {phi_path(y0):.3f}')
    np.save('basins.npy', {'phi': PH, **res}, allow_pickle=True)
    open('basins.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')

    plt.rcParams.update({'font.size': 12, 'axes.labelsize': 14})
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    sty = {'punish': ('tab:red', 'o', 'punishment'),
           'y0=0.5': ('tab:green', 's', 'state dependent rule, switch at $y_0=0.5$'),
           f'y0={round(y0_opt,3)}': ('tab:blue', 'D', f'state dependent rule, switch at $y_0=1/3$'),
           'y0=0.2': ('tab:purple', '^', 'state dependent rule, switch at $y_0=0.2$')}
    for k, v in res.items():
        c, mk, lab = sty[k]
        ax.plot(PH, v, marker=mk, color=c, lw=2.2, label=lab)
    ax.axvline(P['cR']/P['bR'], color='gray', ls=':', lw=1.6)
    ax.axvline((P['cR'] + L)/P['bR'], color='gray', ls='--', lw=1.6)
    ax.axvline(phi_path(y0_opt), color='tab:blue', ls='-.', lw=1.4)
    ax.set_xlabel('$\\phi$'); ax.set_ylabel('share of initial states reaching the target')
    ax.set_ylim(-0.03, 1.03); ax.legend(fontsize=10.5, loc='lower right')
    fig.tight_layout(); fig.savefig('figQ_basins.png', dpi=200); plt.close(fig)
    print('wrote figQ_basins.png basins.txt')
