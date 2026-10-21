"""Pure reward (w = 0): compliance and audit cycle instead of converging.

Checks the statement in SI, Section 4: under pure reward the target state is stable only if
phi (b_R - b_R5) > c_R + r_m + r_h (phi > 0.9 at the illustrative values); for 0.1 < phi < 0.9
the face x = 1 carries a conserved quantity and the orbits are closed. Prints the share of 200
random initial states that reach (1,1,1) for several values of phi, and the variation of the
conserved quantity along one orbit at phi = 0.5.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp

sys.stdout.reconfigure(encoding='utf-8')
bU, th, bP, cP, gm, gh, um, uh = 1, 0.3, 1, 0.3, 0.4, 0.4, 0.8, 0.8
rm, rh, bR, bRf, cR, v, lam = 0.8, 0.8, 2, 0, 0.2, 0.4, 0.4


def f(t, s, phi):
    x, y, z = s
    m = 1 - phi * (1 - z)
    dU = bU * m * (y + th * (1 - y))
    dC = -(cP + gm) - gh * m * x + rm * z + rh * x * z
    dR = -cR + phi * (bR - bRf) * x - y * (rm + rh * x)
    return [x * (1 - x) * dU, y * (1 - y) * dC, z * (1 - z) * dR]


def H(y, z, phi):
    R0, R1 = phi * (bR - bRf) - cR, phi * (bR - bRf) - cR - (rm + rh)
    C0, C1 = -(cP + gm) - gh * (1 - phi), -(cP + gm + gh) + rm + rh
    return R0 * np.log(y) - R1 * np.log(1 - y) - C0 * np.log(z) + C1 * np.log(1 - z)


rng = np.random.default_rng(1)
for phi in (0.05, 0.2, 0.5, 0.85, 0.95):
    den = 0
    for _ in range(200):
        s0 = rng.uniform(0.02, 0.98, 3)
        sol = solve_ivp(f, (0, 3000), s0, args=(phi,), rtol=1e-9, atol=1e-12)
        if np.linalg.norm(sol.y[:, -1] - 1) < 0.05:
            den += 1
    print(f'phi={phi}: {den}/200 initial states reach (1,1,1)')

phi = 0.5
sol = solve_ivp(f, (0, 6000), [1 - 1e-12, 0.3, 0.5], args=(phi,), rtol=1e-11, atol=1e-13,
                dense_output=True)
tt = np.linspace(3000, 6000, 20000)
y, z = sol.sol(tt)[1], sol.sol(tt)[2]
h = H(y, z, phi)
print(f'phi=0.5, x=1: y in [{y.min():.3f},{y.max():.3f}], z in [{z.min():.3f},{z.max():.3f}],'
      f' max deviation of H {np.ptp(h):.2e}')
print('stability threshold of pure reward (cR+rm+rh)/(bR-bRf) =', (cR + rm + rh) / (bR - bRf))
