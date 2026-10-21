"""Reproduces Figures 2 (replicator dynamics) and 3 (finite population heatmaps) of Balabanova
et al. (2025) with media_model.py. Writes fig2_reproduction.png, fig3_reproduction.png and
fig3_data.npy.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import media_model as M

base2 = dict(bI=1.0, cI=0.5, cW=1.0, pW=0.5, bU=4.0, eps=0.2, bP=4.0, cP=0.5, u=1.5,
             bR=4.0, cR=1.0, v=0.5, bfo=1.0)
panels = [(1.0,1.0), (1.0,5.0), (5.0,5.0), (5.0,1.0)]
colors = dict(com='tab:green', usr='tab:red', dev='tab:blue', reg='black')
labels = ['C Commentariat','CT User','C Developer','C Regulator']

fig, axes = plt.subplots(4, 2, figsize=(9, 12))
for r,(bI_,cR_) in enumerate(panels):
    for c,model in enumerate((1,2)):
        p = dict(base2, bI=bI_, cR=cR_)
        sol = solve_ivp(M.rhs, [0,2000], [0.5]*4, args=(model,p), rtol=1e-8, atol=1e-10,
                        dense_output=True)
        t = np.linspace(0,2000,800); Y = sol.sol(t)
        ax = axes[r,c]
        for i,(k,lab) in enumerate(zip(colors, labels)):
            ax.plot(t, Y[i], color=list(colors.values())[i], label=lab)
        ax.set_ylim(-0.02,1.02); ax.set_title(f"Model {'I' if model==1 else 'II'}: $c_I$=0.5, $b_I$={bI_:.0f}, $c_R$={cR_:.0f}", fontsize=9)
        if r==0 and c==0: ax.legend(fontsize=7)
        if r==3: ax.set_xlabel('Iteration')
        if c==0: ax.set_ylabel('Frequency')
fig.suptitle('Reproduction of media-paper Fig. 2 (infinite population, replicator ODE)', y=0.995)
fig.tight_layout()
fig.savefig('fig2_reproduction.png', dpi=130)
print('fig2 done')

base3 = dict(bI=0.0, cI=5.0, cW=0.0, pW=0.5, bU=4.0, eps=0.2, bP=4.0, cP=0.5, u=1.5,
             bR=4.0, cR=0.5, v=0.0, bfo=1.0)
grid = np.linspace(0, 10, 41)
rows = [(1,0.2),(2,0.2),(1,-0.1),(2,-0.1)]
res = np.zeros((4, 4, len(grid), len(grid)))
for ri,(model,eps_) in enumerate(rows):
    for icw, cw in enumerate(grid):
        for ibi, bi in enumerate(grid):
            p = dict(base3, bI=bi, cW=cw, eps=eps_)
            res[ri,:,icw,ibi] = M.coop_freqs(model, p, Z=100, beta=0.1)
    print(f'row {ri+1}/4 done')

np.save('fig3_data.npy', res)
titles = ['Avg. C Commentariat','Avg. CT Users','Avg. C Developers','Avg. C Regulators']
rowlab = ['Model I\n$\\epsilon$=0.2','Model II\n$\\epsilon$=0.2','Model I\n$\\epsilon$=-0.1','Model II\n$\\epsilon$=-0.1']
fig, axes = plt.subplots(4, 4, figsize=(13, 11), constrained_layout=True)
for ri in range(4):
    for pi in range(4):
        ax = axes[ri,pi]
        im = ax.imshow(res[ri,pi], origin='lower', extent=[0,10,0,10], vmin=0, vmax=1,
                       cmap='coolwarm', aspect='auto')
        if ri==0: ax.set_title(titles[pi], fontsize=10)
        if pi==0: ax.set_ylabel(rowlab[ri]+'\n$c_W$', fontsize=9)
        if ri==3: ax.set_xlabel('$b_I$')
fig.colorbar(im, ax=axes, orientation='horizontal', fraction=0.04, pad=0.03, label='Frequency')
fig.suptitle('Reproduction of media-paper Fig. 3 (finite population, Z=100, $\\beta$=0.1)')
fig.savefig('fig3_reproduction.png', dpi=130)
print('fig3 done')
