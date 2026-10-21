"""Symbolic check: the fitness differences derived from payoff tables III and IV coincide with
equations (11a-d) (Model I) and (12a-e) (Model II) of Balabanova et al. (2025).

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import sympy as sp
import media_model as M

x,y,z,w = sp.symbols('x y z w')
bI,cI,cW,pW,bU,eps,bP,cP,u,bR,cR,v,bfo = sp.symbols('b_I c_I c_W p_W b_U epsilon b_P c_P u b_R c_R v b_fo')
p = dict(zip(M.PARAM_KEYS, [bI,cI,cW,pW,bU,eps,bP,cP,u,bR,cR,v,bfo]))

paper = {
 1: [ y*pW*(bI+cW) - cI,
      bU*((x-1)*pW*((z-1)*eps+z) + z),
      (x-1)*y*pW*(2*bP-u*w) + y*bP - cP,
      -(x-1)*y*(z-1)*pW*(v-bfo) - cR ],
 2: [ -cI + y*(bI+cW)*pW,
      -bU*(eps*(-1+z)-z)*(w + (-1+2*w)*(-1+x)*pW),
      -cP + u*w*y + u*w*(-1+x)*y*pW,
      -cR + y*(bfo+bR+v*(-1+z)-bfo*z) + (-1+x)*y*(bfo+2*bR+v*(-1+z)-bfo*z)*pW ],
}
names = ['Commentariat','Users','Developers','Regulators']
for model in (1,2):
    print(f"=== Model {model} ===")
    d = M.deltas(model, p, (x,y,z,w))
    for i in range(4):
        diff = sp.simplify(sp.expand(d[i]) - sp.expand(paper[model][i]))
        print(f"  {names[i]:13s}: delta(Table) - eq_paper = {diff}   {'OK' if diff==0 else '<-- MISMATCH'}")
