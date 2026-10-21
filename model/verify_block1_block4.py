"""Step by step check of the fitness functions and of the stability of the vertices.

Every fitness is built as the full combination of the four profiles of the other two populations
and compared with the closed forms; the Jacobian is evaluated at the eight vertices (eigenvalues
and diagonal structure), and the rest point on the face z = 1 is computed.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import sympy as sp

x, y, z = sp.symbols('x y z', real=True)
bU, bP, cP, gm, gh, um, uh = sp.symbols('b_U b_P c_P g_m g_h u_m u_h', positive=True)
bR, bR5, cR, v, lam = sp.symbols('b_R b_R5 c_R v lambda', positive=True)
theta = sp.Symbol('theta', real=True)
phi = sp.Symbol('phi', nonnegative=True)

Uu = {}
Uc = {}
Ur = {}
Uu[('Use','C','A')] = bU;              Uc[('Use','C','A')] = bP - cP;                       Ur[('Use','C','A')] = bR - cR
Uu[('Use','C','L')] = (1-phi)*bU;      Uc[('Use','C','L')] = (1-phi)*bP - cP;               Ur[('Use','C','L')] = phi*bR5 + (1-phi)*bR
Uu[('Use','D','A')] = theta*bU;        Uc[('Use','D','A')] = bP + gm + gh - um - uh;        Ur[('Use','D','A')] = bR - cR - v - lam
Uu[('Use','D','L')] = (1-phi)*theta*bU;Uc[('Use','D','L')] = (1-phi)*bP + gm + (1-phi)*gh;  Ur[('Use','D','L')] = phi*bR5 + (1-phi)*bR
Uu[('Av','C','A')]  = 0;               Uc[('Av','C','A')]  = -cP;                           Ur[('Av','C','A')]  = bR5 - cR
Uu[('Av','C','L')]  = 0;               Uc[('Av','C','L')]  = -cP;                           Ur[('Av','C','L')]  = bR5
Uu[('Av','D','A')]  = 0;               Uc[('Av','D','A')]  = gm - um;                       Ur[('Av','D','A')]  = bR5 - cR
Uu[('Av','D','L')]  = 0;               Uc[('Av','D','L')]  = gm;                            Ur[('Av','D','L')]  = bR5

pc = {'C': y, 'D': 1-y}
ps = {'A': z, 'L': 1-z}
pu = {'Use': x, 'Av': 1-x}

def fit_user(s):  return sp.expand(sum(pc[c]*ps[st]*Uu[(s,c,st)] for c in 'CD' for st in 'AL'))
def fit_conc(s):  return sp.expand(sum(pu[u]*ps[st]*Uc[(u,s,st)] for u in ('Use','Av') for st in 'AL'))
def fit_state(s): return sp.expand(sum(pu[u]*pc[c]*Ur[(u,c,s)] for u in ('Use','Av') for c in 'CD'))

fU_use, fU_av = fit_user('Use'), fit_user('Av')
fC_c,  fC_d  = fit_conc('C'),  fit_conc('D')
fR_a,  fR_l  = fit_state('A'), fit_state('L')

m = 1 - phi*(1-z)
print("=== BLOCK 1: full fitness (4 cells) ===")
print("f^U_Use  =", sp.factor(fU_use))
print("f^C_Comp =", sp.collect(fC_c, x))
print("f^C_Def  =", sp.collect(fC_d, x))
print("f^R_Aud  =", sp.collect(fR_a, x))
print("f^R_Lax  =", sp.collect(fR_l, x))

r1 = bU*(y + theta*(1-y))
r2c, r2d = bP*x - cP, bP*x + gm + gh*x - um*z - uh*x*z
r3a, r3l = bR*x - cR - (v+lam)*x*(1-y), bR*x
print("\n--- Compare with the stated formulas (phi=0) ---")
for name, full, rep in [("f^U_Use", fU_use, r1), ("f^C_Comp", fC_c, r2c), ("f^C_Def", fC_d, r2d),
                        ("f^R_Aud", fR_a, r3a), ("f^R_Lax", fR_l, r3l)]:
    diff = sp.simplify(full.subs(phi, 0) - rep)
    print(f"{name}: full(phi=0) - stated = {diff}   {'OK (=0)' if diff==0 else '<-- MISMATCH'}")

dU = sp.simplify(fU_use - fU_av)
dC = sp.simplify(fC_c - fC_d)
dR = sp.simplify(fR_a - fR_l)
print("\n=== BLOCK 2: fitness differences (general phi) ===")
rep_dU = bU*m*(y + theta*(1-y))
rep_dC = -(cP+gm) - gh*x*m + um*z + uh*x*z
rep_dR = -cR + phi*(bR-bR5)*x - (v+lam)*x*(1-y)
for name, mine, rep in [("dU", dU, rep_dU), ("dC", dC, rep_dC), ("dR", dR, rep_dR)]:
    d = sp.simplify(mine - rep)
    print(f"{name} - report(13) = {d}   {'OK (=0)' if d==0 else '<-- MISMATCH'}")

print("\n=== BLOCK 4: Jacobian eigenvalues at the 8 vertices (phi=0) ===")
dU0, dC0, dR0 = dU.subs(phi,0), dC.subs(phi,0), dR.subs(phi,0)
xdot = x*(1-x)*dU0; ydot = y*(1-y)*dC0; zdot = z*(1-z)*dR0
J = sp.Matrix([xdot, ydot, zdot]).jacobian([x, y, z])
A = cP + gm; g = gm + gh; L = v + lam
expected = {
 (0,0,0): (bU*theta, -A, -cR),
 (1,0,0): (-bU*theta, -(cP+g), -cR-L),
 (0,1,0): (bU, A, -cR),
 (0,0,1): (bU*theta, um-A, cR),
 (1,1,0): (-bU, cP+g, -cR),
 (1,0,1): (-bU*theta, (um+uh)-(cP+g), cR+L),
 (0,1,1): (bU, A-um, cR),
 (1,1,1): (-bU, (cP+g)-(um+uh), cR),
}
for vtx, exp_eigs in expected.items():
    Jv = J.subs({x:vtx[0], y:vtx[1], z:vtx[2]})
    eigs = [sp.simplify(Jv[i,i]) for i in range(3)]
    ok = all(sp.simplify(e - ex) == 0 for e, ex in zip(eigs, exp_eigs))
    print(f"vertex {vtx}: J diag = {eigs}  vs report {list(exp_eigs)}  {'OK' if ok else '<-- MISMATCH'}")
    offdiag_zero = all(sp.simplify(Jv[i,j])==0 for i in range(3) for j in range(3) if i!=j)
    if not offdiag_zero: print(f"   !!! Jacobian NOT diagonal at {vtx}")

print("\n--- Rest point on the face z=1 ---")
solx = sp.solve(sp.Eq(dC0.subs(z,1), 0), x); soly = sp.solve(sp.Eq(dU0.subs(z,1), 0), y)
print("dC=0 tai z=1  -> x* =", sp.simplify(solx[0]), "   report: (c_P+g_m-u_m)/(u_h-g_h)")
print("dU=0          -> y* =", sp.simplify(soly[0]), "   report: theta/(theta-1)")
print("dR(z=1,phi=0) =", sp.simplify(dR0.subs(z,1)), " (<0 => z=1 cannot be stable; the z eigenvalue -dR is positive)")

print("\n=== Thresholds ===")
dR_good = sp.simplify(dR.subs({x:1, y:1}))
print("dR(1,1,.) =", dR_good, " > 0 <=> phi >", sp.solve(sp.Eq(dR_good,0), phi)[0], "  [phi* = cR/(bR-bR5)] ")
dR_defect = sp.simplify(dR.subs({x:1, y:0}))
print("dR(1,0,.) =", dR_defect, " > 0 <=> phi >", sp.solve(sp.Eq(dR_defect,0), phi)[0], "  [phi** = (cR+v+lam)/(bR-bR5)]")
