"""Fixation probability when the payoff difference depends on the number of mutants.

Under the state dependent rule the instrument weight w depends on the compliance share y of the
concessionaire population itself, so during an invasion of that population the payoff difference
Delta f_j changes with the mutant count j. The general birth-death formula
    rho = [1 + sum_{i=1}^{Z-1} exp(-beta * sum_{j=1}^{i} Delta f_j)]^{-1}
is used, with w evaluated at y = j/Z (Comply mutants) or y = 1 - j/Z (Defect mutants). For users
and the State the payoff difference does not depend on j, and the formula reduces to the geometric
one.

Instrument rules: 'punish' (w = 1), 'reward' (w = 0), 'step' (w = 1[y >= y0]),
'logistic' (w = 1 / (1 + exp(-(y - y0)/s))).

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
import botppp_model as B

STATES = B.STATES


def weight(rule, y, y0=0.5, s=0.05):
    if rule == 'punish': return 1.0
    if rule == 'reward': return 0.0
    if rule == 'step': return 1.0 if y >= y0 - 1e-12 else 0.0
    if rule == 'logistic': return 1.0/(1.0 + np.exp(-(y - y0)/s))
    raise ValueError(rule)


def fit_profile(p, prof, w):
    """Fitness of each of the two strategies of each population at a pure profile of the
    other two populations, for a given weight w. Returns T[(u,c,s)] payoffs."""
    return B.payoff_cells(p, w=w)


def conc_delta(p, x, z, comply_share, mutant, w):
    """Payoff of the mutant strategy minus the resident strategy for a concessionaire,
    with users at pure x and State at pure z."""
    T = B.payoff_cells(p, w=w)
    pc = T[(x, 1, z)][1]
    pd = T[(x, 0, z)][1]
    return (pc - pd) if mutant == 1 else (pd - pc)


def rho_general(dfs, beta):
    """dfs[j-1] = Delta f_j, j = 1..Z-1."""
    cum = np.cumsum(-beta*np.asarray(dfs, float))
    m = cum.max() if cum.size else 0.0
    if m > 700:
        return 0.0
    return 1.0/(1.0 + np.exp(cum).sum())


def transition_matrix(p, Z=50, beta=0.5, rule='punish', y0=0.5, s=0.05, general=True):
    idx = {st: i for i, st in enumerate(STATES)}
    L = np.zeros((8, 8))
    for st in STATES:
        i = idx[st]
        x, y, z = st
        for pop in range(3):
            t = list(st); t[pop] = 1 - t[pop]; t = tuple(t)
            if pop == 1 and general:
                mut = t[1]
                dfs = []
                for j in range(1, Z):
                    yy = j/Z if mut == 1 else 1 - j/Z
                    w = weight(rule, yy, y0, s)
                    dfs.append(conc_delta(p, x, z, yy, mut, w))
                r = rho_general(dfs, beta)
            else:
                w = weight(rule, float(y), y0, s)
                T = B.payoff_cells(p, w=w)
                df = T[t][pop] - T[st][pop]
                r = B.fixation(df, Z, beta)
            L[i, idx[t]] = r/3.0
        L[i, i] = 1.0 - L[i].sum()
    return L


def stationary(p, **kw):
    L = transition_matrix(p, **kw)
    ev, V = np.linalg.eig(L.T)
    k = np.argmin(np.abs(ev - 1.0))
    pi = np.abs(np.real(V[:, k])); return pi/pi.sum()


def pi_good(p, **kw):
    return stationary(p, **kw)[STATES.index((1, 1, 1))]


def liftoff(p, phis, **kw):
    pg = np.array([pi_good(dict(p, phi=ph), **kw) for ph in phis])
    return (phis[np.argmax(pg > 0.5)] if pg.max() > 0.5 else np.nan), pg


if __name__ == '__main__':
    P = dict(B.DEFAULT, bR5=0.0)
    phis = np.linspace(0, 1, 101)
    print('=== old geometric vs general formula, step rule y0=1/2 ===')
    for Z, beta in [(50, 0.5), (50, 10.0), (25, 0.1), (500, 2.0)]:
        old = np.array([B.pi_good(dict(P, scheme='adaptive', phi=ph), Z=Z, beta=beta) for ph in phis])
        new = np.array([pi_good(dict(P, phi=ph), Z=Z, beta=beta, rule='step') for ph in phis])
        print(f'Z={Z} beta={beta}: max |pi_old - pi_new| = {np.abs(old-new).max():.2e}')
    print('\n=== canonical r=u: concessionaire Delta at corners, w=0 vs w=1 ===')
    for x in (0, 1):
        for z in (0, 1):
            d0 = conc_delta(P, x, z, 0, 1, 0.0); d1 = conc_delta(P, x, z, 0, 1, 1.0)
            print(f'x={x} z={z}: Delta(comply-defect) w=0 {d0:+.3f}  w=1 {d1:+.3f}')
