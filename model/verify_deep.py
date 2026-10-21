"""Second, independent numerical verification (no sympy).

 A. Random parameter draws: fitness differences from the 8 payoff cells vs the closed forms.
 B. Full eigenvalues (numpy.linalg.eig, no diagonal assumption) at the 8 vertices.
 C. Inspector reward: Delta_R = -c_R + x(1-y)(b_fo - v - lambda) + phi (b_R - b_R5) x.
 D. ODE runs: (i) the inspector reward gives closed orbits; (ii) the exit valve gives three
    regimes; (iii) the valve ceiling b_R5* = b_R - c_R at phi = 1.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
from scipy.integrate import solve_ivp

rng = np.random.default_rng(0)

def payoff_cells(p, phi, bfo=0.0):
    """Return three dicts (u, c, s) -> payoff, with u, c, s in {0, 1}
    (Use = 1, Comply = 1, Audit = 1)."""
    bU,bP,cP,gm,gh,um,uh,bR,bR5,cR,v,lam,theta = p
    Uu, Uc, Ur = {}, {}, {}
    Uu[(1,1,1)]=bU;             Uc[(1,1,1)]=bP-cP;              Ur[(1,1,1)]=bR-cR
    Uu[(1,1,0)]=(1-phi)*bU;     Uc[(1,1,0)]=(1-phi)*bP-cP;      Ur[(1,1,0)]=phi*bR5+(1-phi)*bR
    Uu[(1,0,1)]=theta*bU;       Uc[(1,0,1)]=bP+gm+gh-um-uh;     Ur[(1,0,1)]=bR-cR-v-lam+bfo
    Uu[(1,0,0)]=(1-phi)*theta*bU; Uc[(1,0,0)]=(1-phi)*bP+gm+(1-phi)*gh; Ur[(1,0,0)]=phi*bR5+(1-phi)*bR
    Uu[(0,1,1)]=0;              Uc[(0,1,1)]=-cP;                Ur[(0,1,1)]=bR5-cR
    Uu[(0,1,0)]=0;              Uc[(0,1,0)]=-cP;                Ur[(0,1,0)]=bR5
    Uu[(0,0,1)]=0;              Uc[(0,0,1)]=gm-um;              Ur[(0,0,1)]=bR5-cR
    Uu[(0,0,0)]=0;              Uc[(0,0,0)]=gm;                 Ur[(0,0,0)]=bR5
    return Uu, Uc, Ur

def deltas_from_cells(xyz, p, phi, bfo=0.0):
    """Fitness differences computed directly from the 8 cells (no closed form)."""
    x,y,z = xyz
    Uu,Uc,Ur = payoff_cells(p, phi, bfo)
    py = {1:y, 0:1-y}; pz = {1:z, 0:1-z}; px = {1:x, 0:1-x}
    fU = {s: sum(py[c]*pz[st]*Uu[(s,c,st)] for c in (0,1) for st in (0,1)) for s in (0,1)}
    fC = {s: sum(px[u]*pz[st]*Uc[(u,s,st)] for u in (0,1) for st in (0,1)) for s in (0,1)}
    fR = {s: sum(px[u]*py[c]*Ur[(u,c,s)]  for u in (0,1) for c in (0,1)) for s in (0,1)}
    return fU[1]-fU[0], fC[1]-fC[0], fR[1]-fR[0]

def deltas_closed(xyz, p, phi, bfo=0.0):
    bU,bP,cP,gm,gh,um,uh,bR,bR5,cR,v,lam,theta = p
    x,y,z = xyz
    m = 1 - phi*(1-z)
    dU = bU*m*(y + theta*(1-y))
    dC = -(cP+gm) - gh*x*m + um*z + uh*x*z
    dR = -cR + phi*(bR-bR5)*x - (v+lam)*x*(1-y) + bfo*x*(1-y)*1.0
    return dU, dC, dR

print("=== A. Numeric cross-check (2000 random draws, general phi, including bfo) ===")
maxerr = 0.0
for _ in range(2000):
    p = rng.uniform(0.05, 2.0, size=13); p[12] = rng.uniform(-1.5, 1.0)
    phi = rng.uniform(0,1); bfo = rng.uniform(0,2)
    xyz = rng.uniform(0,1,3)
    d1 = np.array(deltas_from_cells(xyz, p, phi, bfo))
    d2 = np.array(deltas_closed(xyz, p, phi, bfo))
    maxerr = max(maxerr, np.abs(d1-d2).max())
print(f"max |delta_cells - delta_closed| = {maxerr:.3e}  {'OK' if maxerr < 1e-12 else '<-- MISMATCH'}")

print("\n=== B. Full eigenvalues (numerical Jacobian) at the 8 vertices, 200 draws ===")
def rhs(t, s, p, phi, bfo=0.0):
    x,y,z = s
    dU,dC,dR = deltas_from_cells((x,y,z), p, phi, bfo)
    return [x*(1-x)*dU, y*(1-y)*dC, z*(1-z)*dR]

def num_jac(s, p, phi, bfo=0.0, h=1e-6):
    J = np.zeros((3,3))
    for j in range(3):
        e = np.zeros(3); e[j] = h
        J[:,j] = (np.array(rhs(0, s+e, p, phi, bfo)) - np.array(rhs(0, s-e, p, phi, bfo)))/(2*h)
    return J

def diag_pred(vtx, p, phi, bfo=0.0):
    x,y,z = vtx
    dU,dC,dR = deltas_from_cells(vtx, p, phi, bfo)
    return np.array([(1-2*x)*dU, (1-2*y)*dC, (1-2*z)*dR])

bad = 0
for _ in range(200):
    p = rng.uniform(0.05, 2.0, size=13); p[12] = rng.uniform(-1.5, 1.0)
    phi = rng.uniform(0,1)
    for vtx in [(a,b,c) for a in (0,1) for b in (0,1) for c in (0,1)]:
        J = num_jac(np.array(vtx,float), p, phi)
        eig = np.sort(np.linalg.eigvals(J).real)
        pred = np.sort(diag_pred(vtx, p, phi))
        offd = np.abs(J - np.diag(np.diag(J))).max()
        if np.abs(eig-pred).max() > 1e-4 or offd > 1e-4: bad += 1
print(f"eigenvalue/off-diagonal deviations: {bad} / 1600 cases  {'OK' if bad==0 else '<-- MISMATCH'}")

print("\n=== C/D-i. Reward variant: limit cycle? (bfo > v+lam, theta > 0) ===")
p0 = [1.0, 1.0, 0.3, 0.4, 0.4, 0.8, 0.8, 1.0, 0.0, 0.3, 0.5, 0.4, 0.3]
bfo = 2.5
sol = solve_ivp(rhs, [0, 3000], [0.8,0.8,0.8], args=(p0, 0.0, bfo), dense_output=False,
                rtol=1e-10, atol=1e-12, max_step=1.0)
y_tail = sol.y[1, sol.t > 2000]; z_tail = sol.y[2, sol.t > 2000]
amp_y = y_tail.max()-y_tail.min(); amp_z = z_tail.max()-z_tail.min()
print(f"bfo=2.5: t>2000 amplitude y = {amp_y:.3f}, z = {amp_z:.3f}, end = {np.round(sol.y[:,-1],3)}")
print("   ->", "SUSTAINED OSCILLATION around (y*,z*)=(%.3f,%.3f) OK" % (1-0.3/1.6, 1.1/1.6)
      if (amp_y>0.05 and amp_z>0.05) else "<-- NO oscillation, converges")

sol2 = solve_ivp(rhs, [0, 3000], [0.8,0.8,0.8], args=(p0, 0.0, 0.5), rtol=1e-9, atol=1e-11)
print(f"bfo=0.5 < v+lam: final (x,y,z) = {np.round(sol2.y[:,-1],3)} -> expected z->0, y->0")

print("\n=== D-ii. Exit valve (phi=1): 3 regimes from the same start (0.8,0.8,0.8) ===")
for bR5, label in [(None,"baseline phi=0"), (0.0,"valve OPEN (bR5=0 < 0.7)"), (0.9,"valve CLOSED (bR5=0.9 > 0.7)")]:
    p = list(p0)
    if bR5 is None:
        phi = 0.0
    else:
        p[8] = bR5; phi = 1.0
    s = solve_ivp(rhs, [0, 5000], [0.8,0.8,0.8], args=(p, phi, 0.0), rtol=1e-9, atol=1e-11)
    print(f"{label:28s} -> (x,y,z) = {np.round(s.y[:,-1],4)}")

print("\n=== D-iii. Sweep b_R5 around the ceiling bR-cR=0.7 (phi=1): is the target stable? ===")
for bR5 in [0.6, 0.65, 0.69, 0.71, 0.75, 0.8]:
    p = list(p0); p[8] = bR5
    lam3 = -deltas_from_cells((1,1,1), p, 1.0)[2]
    s = solve_ivp(rhs, [0, 8000], [0.95,0.95,0.95], args=(p, 1.0, 0.0), rtol=1e-9, atol=1e-11)
    end = np.round(s.y[:,-1], 3)
    print(f"bR5={bR5:.2f}: lambda_z(1,1,1) = {lam3:+.3f} ({'stable' if lam3<0 else 'UNSTABLE'}), ODE from (.95,.95,.95) -> {end}")
