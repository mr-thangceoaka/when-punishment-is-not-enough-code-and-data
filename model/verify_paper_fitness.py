"""Independent check of the payoff table as printed in the paper.

Rebuilds every fitness as the full sum over the four profiles of the other two populations,
compares it with the closed forms in the paper (fitness functions, fitness differences, thresholds
phi* and phi**), and cross-checks numerically against botppp_model.fitness.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import sympy as sp
import numpy as np
import botppp_model as B

bU,bP,cP,gm,gh,um,uh = sp.symbols('b_U b_P c_P g_m g_h u_m u_h', positive=True)
bR,bR5,cR,v,lam,theta,phi = sp.symbols('b_R b_R5 c_R v lambda theta phi')
x,y,z = sp.symbols('x y z')

T = {
 (1,1,1): (bU,                bP-cP,                     bR-cR),
 (1,1,0): ((1-phi)*bU,        (1-phi)*bP-cP,             phi*bR5+(1-phi)*bR),
 (1,0,1): (theta*bU,          bP+gm+gh-um-uh,            bR-cR-v-lam),
 (1,0,0): ((1-phi)*theta*bU,  (1-phi)*bP+gm+(1-phi)*gh,  phi*bR5+(1-phi)*bR),
 (0,1,1): (0,                 -cP,                       bR5-cR),
 (0,1,0): (0,                 -cP,                       bR5),
 (0,0,1): (0,                 gm-um,                     bR5-cR),
 (0,0,0): (0,                 gm,                        bR5),
}

pr = {0:{1:x,0:1-x}, 1:{1:y,0:1-y}, 2:{1:z,0:1-z}}
def fit(pop, s):
    others = [j for j in range(3) if j != pop]
    tot = 0
    for a in (0,1):
        for b_ in (0,1):
            prof = [None]*3; prof[pop]=s
            prof[others[0]], prof[others[1]] = a, b_
            tot += pr[others[0]][a]*pr[others[1]][b_]*T[tuple(prof)][pop]
    return sp.expand(tot)

m = 1 - phi*(1-z)

paper = {
 'fU_Use'   : bU*m*(y+theta*(1-y)),
 'fU_Avoid' : sp.Integer(0),
 'fC_Comply': bP*m*x - cP,
 'fC_Defect': bP*m*x + gm + gh*m*x - um*z - uh*x*z,
 'fR_Audit' : x*bR + (1-x)*bR5 - cR - (v+lam)*x*(1-y),
 'fR_Lax'   : x*((1-phi)*bR + phi*bR5) + (1-x)*bR5,
}
calc = {
 'fU_Use': fit(0,1), 'fU_Avoid': fit(0,0),
 'fC_Comply': fit(1,1), 'fC_Defect': fit(1,0),
 'fR_Audit': fit(2,1), 'fR_Lax': fit(2,0),
}
print('=== 1) Paper payoff table -> fitness vs closed forms in the paper ===')
ok = True
for k in paper:
    d = sp.simplify(calc[k] - paper[k])
    print(f'  {k:11s}: hieu = {d}  {"OK" if d==0 else "WRONG!!!"}')
    ok &= (d == 0)

print('=== 2) Differences (eqs dU, dC, dR in the paper) ===')
dU_p = bU*m*(y+theta*(1-y))
dC_p = -(cP+gm) - gh*m*x + um*z + uh*x*z
dR_p = -cR + phi*(bR-bR5)*x - (v+lam)*x*(1-y)
for name, lhs, rhs in [('dU', calc['fU_Use']-calc['fU_Avoid'], dU_p),
                       ('dC', calc['fC_Comply']-calc['fC_Defect'], dC_p),
                       ('dR', calc['fR_Audit']-calc['fR_Lax'], dR_p)]:
    d = sp.simplify(lhs - rhs)
    print(f'  {name}: hieu = {d}  {"OK" if d==0 else "WRONG!!!"}')
    ok &= (d == 0)

print('=== 3) Thresholds phi*, phi** (eqs in the paper) ===')
dR = calc['fR_Audit']-calc['fR_Lax']
t1 = sp.solve(sp.Eq(dR.subs({x:1,y:1}), 0), phi)[0]
t2 = sp.solve(sp.Eq(dR.subs({x:1,y:0}), 0), phi)[0]
e1 = sp.simplify(t1 - cR/(bR-bR5)); e2 = sp.simplify(t2 - (cR+v+lam)/(bR-bR5))
print(f'  phi*  - cR/(bR-bR5)          = {e1}  {"OK" if e1==0 else "WRONG!!!"}')
print(f'  phi** - (cR+v+lam)/(bR-bR5)  = {e2}  {"OK" if e2==0 else "WRONG!!!"}')
ok &= (e1 == 0 and e2 == 0)

print('=== 4) Numerical cross-check with botppp_model.fitness (simulation code) ===')
rng = np.random.default_rng(7)
syms = (bU,bP,cP,gm,gh,um,uh,bR,bR5,cR,v,lam,theta,phi,x,y,z)
fns = {k: sp.lambdify(syms, calc[k], 'numpy') for k in calc}
maxdev = 0.0
for _ in range(500):
    vals = dict(bU=rng.uniform(.1,2), bP=rng.uniform(.1,2), cP=rng.uniform(.05,1),
                gm=rng.uniform(.05,1), gh=rng.uniform(.05,1), um=rng.uniform(.05,2),
                uh=rng.uniform(.05,2), bR=rng.uniform(.5,3), bR5=rng.uniform(0,1.5),
                cR=rng.uniform(.05,1), v=rng.uniform(.05,1), lam=rng.uniform(.05,1),
                theta=rng.uniform(-1,1), phi=rng.uniform(0,1))
    fx,fy,fz = rng.uniform(0,1,3)
    F = B.fitness(vals, (fx,fy,fz))
    args = tuple(vals[k] for k in ('bU','bP','cP','gm','gh','um','uh','bR','bR5',
                                   'cR','v','lam','theta','phi')) + (fx,fy,fz)
    for (pop,s), key in [((0,1),'fU_Use'),((0,0),'fU_Avoid'),((1,1),'fC_Comply'),
                         ((1,0),'fC_Defect'),((2,1),'fR_Audit'),((2,0),'fR_Lax')]:
        maxdev = max(maxdev, abs(F[pop][s] - fns[key](*args)))
print(f'  500 draws x 6 fitness: max |paper - code| = {maxdev:.2e}  '
      f'{"OK" if maxdev < 1e-12 else "WRONG!!!"}')
ok &= (maxdev < 1e-12)

print()
print('==> ALL MATCH.' if ok else '==> MISMATCH, SEE THE LOG ABOVE!')
