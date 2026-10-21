"""Extensions of the model (SI, Sections 5 and 6).

(A) Users who observe the operator's performance: a fraction phi_q of users leave when the operator
    defects, whatever the oversight, in addition to the fraction phi (phi_o) who leave when
    oversight is lax. With phi_q = 0 the payoffs equal botppp_model.payoff_cells.
(B) Social welfare per period W, counting only real gains and costs: tolls, fines, illicit gains
    and rewards are transfers between the parties and do not enter W.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
import botppp_model as B

EXT_DEFAULT = dict(phiq=0.0,
                   Ds=2.0,
                   dq=0.7,
                   kc=0.5,
                   mu=0.3,
                   eta=0.5)


def payoff_cells_ext(p, w=1.0):
    q = {**B.DEFAULT, **EXT_DEFAULT, **p}
    bU, bP, cP, gm, gh, um, uh = q['bU'], q['bP'], q['cP'], q['gm'], q['gh'], q['um'], q['uh']
    bR, bR5, cR, v, lam, th, phi, bfo = q['bR'], q['bR5'], q['cR'], q['v'], q['lam'], q['theta'], q['phi'], q['bfo']
    rm, rh, pq = q['rm'], q['rh'], q['phiq']
    assert phi + pq <= 1 + 1e-12
    T = {}
    T[(1, 1, 1)] = (bU, bP - cP + (1 - w)*(rm + rh), bR - cR - (1 - w)*(rm + rh))
    T[(1, 1, 0)] = ((1 - phi)*bU, (1 - phi)*bP - cP, phi*bR5 + (1 - phi)*bR)
    T[(1, 0, 1)] = ((1 - pq)*th*bU, (1 - pq)*(bP + gh - w*uh) + gm - w*um,
                    pq*bR5 + (1 - pq)*bR - cR - w*(v + lam) + bfo)
    T[(1, 0, 0)] = ((1 - phi - pq)*th*bU, (1 - phi - pq)*(bP + gh) + gm,
                    (phi + pq)*bR5 + (1 - phi - pq)*bR)
    T[(0, 1, 1)] = (0.0, -cP + (1 - w)*rm, bR5 - cR - (1 - w)*rm)
    T[(0, 1, 0)] = (0.0, -cP, bR5)
    T[(0, 0, 1)] = (0.0, gm - w*um, bR5 - cR)
    T[(0, 0, 0)] = (0.0, gm, bR5)
    return T


def gradients(p, x, y, z, w=1.0):
    """Fitness differences (Use-Avoid, Comply-Defect, Audit-Lax) from the four cell sums."""
    T = payoff_cells_ext(p, w)
    pr = [{1: x, 0: 1 - x}, {1: y, 0: 1 - y}, {1: z, 0: 1 - z}]
    d = []
    for i in range(3):
        o = [j for j in range(3) if j != i]
        f = {}
        for s in (0, 1):
            tot = 0.0
            for a in (0, 1):
                for b in (0, 1):
                    prof = [None]*3; prof[i] = s; prof[o[0]], prof[o[1]] = a, b
                    tot += pr[o[0]][a]*pr[o[1]][b]*T[tuple(prof)][i]
            f[s] = tot
        d.append(f[1] - f[0])
    return np.array(d)


def welfare(p, x, y, z, w=1.0):
    """Social welfare W(x, y, z) and the reward outlay per period."""
    q = {**B.DEFAULT, **EXT_DEFAULT, **p}
    phi, pq = q['phi'], q['phiq']
    mo = 1 - phi*(1 - z)
    TC = x*y*mo
    TD = x*(1 - y)*(mo - pq)
    sA = 1 - TC - TD
    rew = z*(1 - w)*y*(q['rm'] + q['rh']*x)
    W = (q['Ds']*TC + (q['Ds'] - q['dq'])*TD - q['kc']*sA**2 - q['cP']*y
         - z*q['cR'] - z*x*(1 - y)*w*(q['v'] + q['eta']*q['lam']) - q['mu']*rew)
    return W, rew


if __name__ == '__main__':
    rng = np.random.default_rng(3)
    mx = 0.0
    for _ in range(500):
        p = dict(phi=rng.random(), bR5=rng.random()*1.5, bfo=rng.random())
        w = rng.random()
        A = payoff_cells_ext(dict(p, phiq=0.0), w); Bt = B.payoff_cells(p, w)
        mx = max(mx, max(abs(a - b) for k in A for a, b in zip(A[k], Bt[k])))
    print(f'phiq=0 reproduces the main payoff table: max deviation {mx:.1e}')
    mx = 0.0
    for _ in range(500):
        ph = rng.random()*0.6; pq = rng.random()*(1 - ph); x, y, z = rng.random(3)
        p = dict(B.DEFAULT, phi=ph, phiq=pq, bR5=rng.random())
        dU, dC, dR = gradients(p, x, y, z)
        mo = 1 - ph*(1 - z)
        dU_c = p['bU']*(y*mo + p['theta']*(1 - y)*(mo - pq))
        dC_c = (-(p['cP'] + p['gm']) + x*p['bP']*pq - x*p['gh']*(mo - pq)
                + p['um']*z + p['uh']*x*z*(1 - pq))
        dR_c = -p['cR'] + ph*(p['bR'] - p['bR5'])*x - (p['v'] + p['lam'])*x*(1 - y)
        mx = max(mx, abs(dU - dU_c), abs(dC - dC_c), abs(dR - dR_c))
    print(f'closed form gradients with phiq: max deviation {mx:.1e}')
