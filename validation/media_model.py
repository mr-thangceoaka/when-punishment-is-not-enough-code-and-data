"""Reimplementation of the model of Balabanova et al. (2025), "Media and responsible AI
governance: a game-theoretic and LLM analysis", arXiv:2503.09858.

Four populations: x = commentariat (C/D), y = users (CT/N), z = developers (C/D), w = regulators
(C/D); 1 = cooperate or conditional trust, 0 = defect or never. Payoff tables III (Model I) and IV
(Model II), fitness as the full sum over the 8 profiles of the other three populations, replicator
dynamics, and the finite population Markov chain on the 16 monomorphic states. Used to validate the
pipeline before applying it to the concession game.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np

PARAM_KEYS = ['bI','cI','cW','pW','bU','eps','bP','cP','u','bR','cR','v','bfo']

def payoff_table(model, p):
    """Return a dict (sc, su, sd, sr) -> (P_com, P_user, P_dev, P_reg)."""
    bI,cI,cW,pW,bU,eps,bP,cP,u,bR,cR,v,bfo = (p[k] for k in PARAM_KEYS)
    comD = (1-pW)*bI - pW*cW
    T = {}
    if model == 1:
        T[(1,1,1,1)] = (bI-cI, bU, bP-cP, bR-cR)
        T[(1,1,1,0)] = (bI-cI, bU, bP-cP, bR)
        T[(1,1,0,1)] = (bI-cI, 0.0, 0.0, -cR)
        T[(1,1,0,0)] = (bI-cI, 0.0, 0.0, 0.0)
        T[(0,1,1,1)] = (comD, (1-pW)*bU, (1-pW)*bP-cP, (1-pW)*bR-cR)
        T[(0,1,1,0)] = (comD, (1-pW)*bU, (1-pW)*bP-cP, (1-pW)*bR)
        T[(0,1,0,1)] = (comD, pW*eps*bU, pW*(bP-u), pW*(bR+bfo-v)-cR)
        T[(0,1,0,0)] = (comD, pW*eps*bU, pW*bP, pW*bR)
    else:
        T[(1,1,1,1)] = (bI-cI, bU, bP-cP, bR-cR)
        T[(1,1,1,0)] = (bI-cI, 0.0, -cP, 0.0)
        T[(1,1,0,1)] = (bI-cI, eps*bU, bP-u, bR-cR-v+bfo)
        T[(1,1,0,0)] = (bI-cI, 0.0, 0.0, 0.0)
        T[(0,1,1,1)] = (comD, (1-pW)*bU, (1-pW)*bP-cP, (1-pW)*bR-cR)
        T[(0,1,1,0)] = (comD, pW*bU, pW*bP-cP, pW*bR)
        T[(0,1,0,1)] = (comD, (1-pW)*eps*bU, (1-pW)*(bP-u), (1-pW)*(bR-cR+bfo-v)-pW*cR)
        T[(0,1,0,0)] = (comD, pW*eps*bU, pW*bP, pW*bR)
    for sc in (0,1):
        for sd in (0,1):
            for sr in (0,1):
                pcom = -cI if sc==1 else 0.0
                pdev = -cP if sd==1 else 0.0
                preg = -cR if sr==1 else 0.0
                T[(sc,0,sd,sr)] = (pcom, 0.0, pdev, preg)
    return T

def fitness(model, p, freqs):
    """freqs = (x, y, z, w). Return fit[i][s]: fitness of strategy s in
    population i, the full sum over the 8 strategy profiles of the three other populations
    (eqs 1-4 of the paper)."""
    T = payoff_table(model, p)
    x,y,z,w = freqs
    pr = [ {1:x,0:1-x}, {1:y,0:1-y}, {1:z,0:1-z}, {1:w,0:1-w} ]
    fit = [dict(), dict(), dict(), dict()]
    for i in range(4):
        others = [j for j in range(4) if j != i]
        for s in (0,1):
            tot = 0.0
            for a in (0,1):
                for b in (0,1):
                    for c in (0,1):
                        prof = [None]*4
                        prof[i] = s
                        prof[others[0]], prof[others[1]], prof[others[2]] = a,b,c
                        tot += pr[others[0]][a]*pr[others[1]][b]*pr[others[2]][c] * T[tuple(prof)][i]
            fit[i][s] = tot
    return fit

def deltas(model, p, freqs):
    f = fitness(model, p, freqs)
    return np.array([f[i][1]-f[i][0] for i in range(4)])

def rhs(t, s, model, p):
    d = deltas(model, p, s)
    return [si*(1-si)*di for si,di in zip(s,d)]

def fixation(delta_f, Z, beta):
    """Fixation probability of a single mutant, delta_f = f_mut - f_res (constant
    in k, because the 4-player game takes one individual from each other population)."""
    a = beta*delta_f
    if abs(a) < 1e-12:
        return 1.0/Z
    with np.errstate(over='ignore'):
        num = -np.expm1(-a)
        den = -np.expm1(-Z*a)
    if not np.isfinite(den) or den == 0:
        return max(num, 0.0) if a > 0 else 0.0
    return num/den

def stationary(model, p, Z=100, beta=0.1):
    """Markov chain on the 16 monomorphic states (eqs 5-7). Return pi (16,)
    and the list of states."""
    states = [(a,b,c,d) for a in (0,1) for b in (0,1) for c in (0,1) for d in (0,1)]
    idx = {s:i for i,s in enumerate(states)}
    n = len(states)
    L = np.zeros((n,n))
    for s in states:
        i = idx[s]
        for pop in range(4):
            t = list(s); t[pop] = 1-t[pop]; t = tuple(t)
            freqs = [float(v) for v in s]
            f = fitness(model, p, freqs)
            df = f[pop][t[pop]] - f[pop][s[pop]]
            L[i, idx[t]] = fixation(df, Z, beta)/4.0
        L[i,i] = 1.0 - L[i].sum()
    evals, evecs = np.linalg.eig(L.T)
    k = np.argmin(np.abs(evals-1.0))
    pi = np.real(evecs[:,k]); pi = np.abs(pi); pi /= pi.sum()
    return pi, states

def coop_freqs(model, p, Z=100, beta=0.1):
    """Long run frequency of cooperation in each population: the sum of pi
    over the states with strategy C."""
    pi, states = stationary(model, p, Z, beta)
    out = np.zeros(4)
    for w_, s in zip(pi, states):
        for i in range(4):
            out[i] += w_*s[i]
    return out
