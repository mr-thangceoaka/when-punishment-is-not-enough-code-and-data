"""Core model: three populations in a highway concession.

x: share of users who take the tolled facility (Use; otherwise Avoid)
y: share of concessionaires who meet the contract and report revenue honestly (Comply; otherwise Defect)
z: share of State agencies that audit (Audit; otherwise Lax)

Defines the illustrative parameter values (DEFAULT), the payoffs of the eight joint profiles
(Table 1 of the paper), the fitness of each strategy as the full sum over the profiles of the other
two populations, the fitness differences and the replicator dynamics, the instrument weight w(y) of
the state dependent rule, and the finite population Markov chain in the limit of rare mutations
(fixation probability, stationary distribution, long run share of the target state).

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np

DEFAULT = dict(bU=1.0, bP=1.0, cP=0.3, gm=0.4, gh=0.4, um=0.8, uh=0.8,
               bR=2.0, bR5=0.0, cR=0.2, v=0.4, lam=0.4, theta=0.3,
               phi=0.0, bfo=0.0,
               scheme='punish', rm=0.8, rh=0.8)

def payoff_cells(p, w=1.0):
    """(u, c, s) -> (U_user, U_conc, U_state); 1 = Use/Comply/Audit.
    w = weight of the instrument under Audit: 1 = pure punishment (default), 0 = pure reward.
    Under reward a complying firm receives r_m (+ r_h if the user is on the facility); the State
    pays that amount (a transfer) and bears no enforcement cost v + lam."""
    q = {**DEFAULT, **p}
    bU,bP,cP,gm,gh,um,uh = q['bU'],q['bP'],q['cP'],q['gm'],q['gh'],q['um'],q['uh']
    bR,bR5,cR,v,lam,theta,phi,bfo = q['bR'],q['bR5'],q['cR'],q['v'],q['lam'],q['theta'],q['phi'],q['bfo']
    rm,rh = q['rm'],q['rh']
    T = {}
    T[(1,1,1)] = (bU,               bP-cP + (1-w)*(rm+rh),  bR-cR - (1-w)*(rm+rh))
    T[(1,1,0)] = ((1-phi)*bU,       (1-phi)*bP-cP,          phi*bR5+(1-phi)*bR)
    T[(1,0,1)] = (theta*bU,         bP+gm+gh-w*(um+uh),     bR-cR-w*(v+lam)+bfo)
    T[(1,0,0)] = ((1-phi)*theta*bU, (1-phi)*bP+gm+(1-phi)*gh, phi*bR5+(1-phi)*bR)
    T[(0,1,1)] = (0.0, -cP + (1-w)*rm,  bR5-cR - (1-w)*rm)
    T[(0,1,0)] = (0.0, -cP,             bR5)
    T[(0,0,1)] = (0.0, gm-w*um,         bR5-cR)
    T[(0,0,0)] = (0.0, gm,              bR5)
    return T

def instrument_weight(q, y, sharp=False):
    """Weight of punishment under each scheme. adaptive: reward when y is low,
    punish when y is high. sharp=True: exact step (used for pure states)."""
    sch = q.get('scheme', 'punish')
    if sch == 'punish': return 1.0
    if sch == 'reward': return 0.0
    if sharp: return 1.0 if y >= 0.5 else 0.0
    return 1.0/(1.0 + np.exp(-(y-0.5)/0.05))

def fitness(p, freqs, sharp=False):
    """fit[i][s] = fitness of strategy s in population i: the full sum over
    the four strategy profiles of the two other populations."""
    q = {**DEFAULT, **p}
    x,y,z = freqs
    T = payoff_cells(p, w=instrument_weight(q, y, sharp))
    pr = [ {1:x,0:1-x}, {1:y,0:1-y}, {1:z,0:1-z} ]
    fit = [dict(), dict(), dict()]
    for i in range(3):
        others = [j for j in range(3) if j != i]
        for s in (0,1):
            tot = 0.0
            for a in (0,1):
                for b in (0,1):
                    prof = [None]*3
                    prof[i] = s
                    prof[others[0]], prof[others[1]] = a,b
                    tot += pr[others[0]][a]*pr[others[1]][b]*T[tuple(prof)][i]
            fit[i][s] = tot
    return fit

def deltas(p, freqs):
    f = fitness(p, freqs)
    return np.array([f[i][1]-f[i][0] for i in range(3)])

def rhs(t, s, p):
    d = deltas(p, s)
    return [si*(1-si)*di for si,di in zip(s,d)]

def deltas_vec(X, Y, Z, q):
    m = 1 - q['phi']*(1-Z)
    dU = q['bU']*m*(Y + q['theta']*(1-Y))
    dC = -(q['cP']+q['gm']) - q['gh']*X*m + q['um']*Z + q['uh']*X*Z
    dR = -q['cR'] + q['phi']*(q['bR']-q['bR5'])*X - (q['v']+q['lam'])*X*(1-Y) + q['bfo']*X*(1-Y)
    return dU, dC, dR

def rk4_batch(X0, q, T=2000.0, dt=0.1):
    """Integrate the replicator dynamics for a batch of initial conditions.
    X0 has shape (N, 3); q is a dict whose values may be scalars or arrays of shape (N,)."""
    S = np.array(X0, dtype=float).copy()
    def f(S):
        X, Y, Z = S[:,0], S[:,1], S[:,2]
        dU, dC, dR = deltas_vec(X, Y, Z, q)
        return np.stack([X*(1-X)*dU, Y*(1-Y)*dC, Z*(1-Z)*dR], axis=1)
    n = int(T/dt)
    for _ in range(n):
        k1 = f(S); k2 = f(S+dt/2*k1); k3 = f(S+dt/2*k2); k4 = f(S+dt*k3)
        S += dt/6*(k1+2*k2+2*k3+k4)
        np.clip(S, 0.0, 1.0, out=S)
    return S

def fixation(delta_f, Z, beta):
    a = beta*delta_f
    if abs(a) < 1e-12: return 1.0/Z
    with np.errstate(over='ignore'):
        num = -np.expm1(-a); den = -np.expm1(-Z*a)
    if not np.isfinite(den) or den == 0:
        return max(num, 0.0) if a > 0 else 0.0
    return num/den

STATES = [(a,b,c) for a in (0,1) for b in (0,1) for c in (0,1)]

def stationary(p, Z=50, beta=0.5):
    idx = {s:i for i,s in enumerate(STATES)}
    L = np.zeros((8,8))
    for s in STATES:
        i = idx[s]
        for pop in range(3):
            t = list(s); t[pop] = 1-t[pop]; t = tuple(t)
            f = fitness(p, [float(v) for v in s], sharp=True)
            df = f[pop][t[pop]] - f[pop][s[pop]]
            L[i, idx[t]] = fixation(df, Z, beta)/3.0
        L[i,i] = 1.0 - L[i].sum()
    evals, evecs = np.linalg.eig(L.T)
    k = np.argmin(np.abs(evals-1.0))
    pi = np.abs(np.real(evecs[:,k])); pi /= pi.sum()
    return pi

def pi_good(p, Z=50, beta=0.5):
    return stationary(p, Z, beta)[STATES.index((1,1,1))]
