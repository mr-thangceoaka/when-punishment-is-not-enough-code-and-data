"""Symbolic checks (sympy) of the formulas stated in the paper and the SI: fitness sums from the
payoff table, fitness differences, the thresholds phi*, phi**, phi_s and phi_dagger, the
extensions and the welfare function. Prints one line per identity and the number that hold.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import sympy as sp

x, y, z, phi, pq = sp.symbols('x y z phi phi_q', real=True)
bU, bP, cP, gm, gh, um, uh, bR, bR5, cR, v, lam, th, rm, rh, w, bfo = sp.symbols(
    'b_U b_P c_P g_m g_h u_m u_h b_R b_R5 c_R v lambda theta r_m r_h w b_fo', real=True)
ok = []
def check(name, expr):
    r = sp.simplify(expr)
    ok.append((name, r == 0)); print(f'{"OK " if r == 0 else "FAIL"} {name}' + ('' if r == 0 else f'   residual: {r}'))

T = {
 (1,1,1): (bU, bP - cP + (1-w)*(rm+rh), bR - cR - (1-w)*(rm+rh)),
 (1,1,0): ((1-phi)*bU, (1-phi)*bP - cP, phi*bR5 + (1-phi)*bR),
 (1,0,1): ((1-pq)*th*bU, (1-pq)*(bP+gh-w*uh) + gm - w*um, pq*bR5 + (1-pq)*bR - cR - w*(v+lam)),
 (1,0,0): ((1-phi-pq)*th*bU, (1-phi-pq)*(bP+gh) + gm, (phi+pq)*bR5 + (1-phi-pq)*bR),
 (0,1,1): (0, -cP + (1-w)*rm, bR5 - cR - (1-w)*rm),
 (0,1,0): (0, -cP, bR5),
 (0,0,1): (0, gm - w*um, bR5 - cR),
 (0,0,0): (0, gm, bR5)}
pr = [{1: x, 0: 1-x}, {1: y, 0: 1-y}, {1: z, 0: 1-z}]
def fit(i, s, TT):
    o = [j for j in range(3) if j != i]; tot = 0
    for a in (0, 1):
        for b in (0, 1):
            prof = [None]*3; prof[i] = s; prof[o[0]], prof[o[1]] = a, b
            tot += pr[o[0]][a]*pr[o[1]][b]*TT[tuple(prof)][i]
    return sp.expand(tot)
def grads(TT):
    return [fit(i, 1, TT) - fit(i, 0, TT) for i in range(3)]

base = {k: tuple(sp.sympify(c).subs({pq: 0, w: 1}) for c in val) for k, val in T.items()}
dU, dC, dR = grads(base)
m = 1 - phi*(1-z)
check('eq (6) Delta_U', dU - bU*m*(y + th*(1-y)))
check('eq (7) Delta_C', dC - (-(cP+gm) - gh*m*x + um*z + uh*x*z))
check('eq (8) Delta_R', dR - (-cR + phi*(bR-bR5)*x - (v+lam)*x*(1-y)))
check('SI eq (S1) baseline Delta_R at phi=0', dR.subs(phi, 0) - (-cR - (v+lam)*x*(1-y)))
check('dDelta_C/dx = u_h z - g_h m', sp.diff(dC, x) - (uh*z - gh*m))
check('SI f_Use', fit(0,1,base) - bU*m*(y+th*(1-y)))
check('SI f_Comply', fit(1,1,base) - (bP*m*x - cP))
check('SI f_Defect', fit(1,0,base) - (bP*m*x + gm + gh*m*x - um*z - uh*x*z))
check('SI f_Audit', fit(2,1,base) - (x*bR + (1-x)*bR5 - cR - (v+lam)*x*(1-y)))
check('SI f_Lax', fit(2,0,base) - (x*((1-phi)*bR + phi*bR5) + (1-x)*bR5))
A, g, L = cP+gm, gm+gh, v+lam
table = {(0,0,0): (bU*th, -A, -cR), (1,0,0): (-bU*th, -(cP+g), -cR-L), (0,1,0): (bU, A, -cR),
         (0,0,1): (bU*th, um-A, cR), (1,1,0): (-bU, cP+g, -cR), (1,0,1): (-bU*th, um+uh-(cP+g), cR+L),
         (0,1,1): (bU, A-um, cR), (1,1,1): (-bU, cP+g-um-uh, cR)}
for vtx, ev in table.items():
    sub = {x: vtx[0], y: vtx[1], z: vtx[2], phi: 0}
    calc = [(1-2*vtx[0])*dU.subs(sub), (1-2*vtx[1])*dC.subs(sub), (1-2*vtx[2])*dR.subs(sub)]
    for k in range(3):
        check(f'SI Table S4 vertex {vtx} lambda_{k+1}', calc[k] - ev[k])
check('Prop 2 Delta_R at (1,1,1)', dR.subs({x:1, y:1}) - (phi*(bR-bR5) - cR))
check('Prop 2 Delta_R at (1,0)', dR.subs({x:1, y:0}) - (phi*(bR-bR5) - cR - (v+lam)))
xs = ((cP+gm) - um)/(uh-gh); ys = th/(th-1)
check('Lemma S1 x*', dC.subs({z:1, x: xs}))
check('Lemma S1 y*', dU.subs({z:1, y: ys}))
dRr = dR + bfo*x*(1-y)
Acoef = -(cP+gm) - gh*(1-phi); Bcoef = um + uh - gh*phi; D = bfo - v - lam; E = D + phi*(bR-bR5) - cR
check('Prop 3 Delta_C on x=1 = A + B z', dC.subs(x,1) - (Acoef + Bcoef*z))
check('Prop 3 Delta_R on x=1 = E - D y', dRr.subs(x,1) - (E - D*y))
H = E*sp.log(y) - (E-D)*sp.log(1-y) - Acoef*sp.log(z) + (Acoef+Bcoef)*sp.log(1-z)
ydot = y*(1-y)*(Acoef+Bcoef*z); zdot = z*(1-z)*(E-D*y)
check('Prop 3 dH/dt = 0', sp.diff(H, y)*ydot + sp.diff(H, z)*zdot)
check('eq bfo threshold at (1,0), phi=0', dRr.subs({x:1, y:0, phi:0}) - (bfo - v - lam - cR))
rule = {k: tuple(sp.sympify(c).subs({pq: 0}) for c in val) for k, val in T.items()}
dUw, dCw, dRw = grads(rule)
check('eq (14) Delta_C under w', dCw - (-(cP+gm) - gh*m*x + w*(um*z+uh*x*z) + (1-w)*(rm*z+rh*x*z)))
check('eq (15) Delta_R under w', dRw - (-cR + phi*(bR-bR5)*x - w*(v+lam)*x*(1-y) - (1-w)*y*(rm+rh*x)))
check('SI eq (S6) c(y) on x=1', dRw.subs(x,1) - (phi*(bR-bR5) - cR - (w*(v+lam)*(1-y) + (1-w)*(rm+rh)*y)))
check('Prop 4 at (1,1) w=1', dRw.subs({x:1, y:1, w:1}) - (phi*(bR-bR5) - cR))
check('Prop 4 at (1,0) w=0', dRw.subs({x:1, y:0, w:0}) - (phi*(bR-bR5) - cR))
check('Prop 5 Delta_C at x=1,z=1', dCw.subs({x:1, z:1}) - (-(cP+gm+gh) + w*(um+uh) + (1-w)*(rm+rh)))
R_ = rm + rh
y0s = L/(L+R_)
check('Prop 5 best y0 equalises the two ends', y0s*R_ - (1-y0s)*L)
check('Prop 5 value at best y0', y0s*R_ - L*R_/(L+R_))
dUq, dCq, dRq = grads({k: tuple(sp.sympify(c).subs({w: 1}) for c in val) for k, val in T.items()})
check('SI eq (S8) Delta_U with phi_q', dUq - bU*(y*m + th*(1-y)*(m - pq)))
check('SI eq (S9) Delta_C with phi_q', dCq - (-(cP+gm) + bP*pq*x - gh*(m-pq)*x + um*z + uh*(1-pq)*x*z))
check('SI eq (S10) Delta_R with phi_q (no phi_q)', dRq - (-cR + phi*(bR-bR5)*x - (v+lam)*x*(1-y)))
check('Prop 6 (ii) Delta_C at x=1,z=0', dCq.subs({x:1, z:0}) - (-(cP+gm) + (bP+gh)*pq - gh*(1-phi)))
Ds, dq, k, mu, eta = sp.symbols('D_s d_q k mu eta', real=True)
def W(xx, yy, zz):
    mo = 1 - phi*(1-zz); TC = xx*yy*mo; TD = xx*(1-yy)*(mo-pq); sA = 1 - TC - TD
    return Ds*TC + (Ds-dq)*TD - k*sA**2 - cP*yy - zz*cR - zz*xx*(1-yy)*(v+eta*lam)
check('Prop 7 W(1,1,1)', W(1,1,1) - (Ds - cP - cR))
check('Prop 7 W(1,1,0)', W(1,1,0) - (Ds*(1-phi) - k*phi**2 - cP))
check('Prop 7 W(1,0,0)', W(1,0,0) - ((Ds-dq)*(1-phi-pq) - k*(phi+pq)**2))
check('Prop 7 W(1,0,1)', W(1,0,1) - ((Ds-dq)*(1-pq) - k*pq**2 - cR - v - eta*lam))
check('Prop 7 W(0,y,z)', W(0,y,z) - (-k - cP*y - cR*z))
check('Prop 7 (b)', (W(1,1,1) - W(1,1,0)) - (phi*(Ds + k*phi) - cR))
check('Prop 7 (c)', (W(1,1,1) - W(1,0,0)) - (Ds*(phi+pq) + dq*(1-phi-pq) + k*(phi+pq)**2 - cP - cR))
check('Prop 7 (d)', (W(1,1,1) - W(1,0,1)) - (Ds*pq + dq*(1-pq) + k*pq**2 + v + eta*lam - cP))
check('Prop 7 (a)', (W(1,1,1) - W(0,0,0)) - (Ds + k - cP - cR))
phis = (2*cR + v + lam)/(2*(bR-bR5))
check('eq (18) phi_s equalises the two State exits', (phis*(bR-bR5) - cR) - (cR + v + lam - phis*(bR-bR5)))
kap = sp.symbols('kappa')
check('phi**_kappa', (dR + kap*(um+uh)*x*(1-y)).subs({x:1, y:0}) - (phi*(bR-bR5) - (cR+v+lam-kap*(um+uh))))
print(f'\n{sum(r for _, r in ok)} of {len(ok)} checks passed')
