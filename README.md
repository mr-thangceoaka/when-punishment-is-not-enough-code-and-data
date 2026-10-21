# When Punishment Is Not Enough: code and data

Code and data for the paper

> Toan Thang Nguyen, The Anh Han, Matjaž Perc and Doanh Nguyen Ngoc,
> **When Punishment Is Not Enough: Incentive Design for Compliance in Highway Public-Private Partnerships** (2026, preprint).

The paper studies an evolutionary game with three populations: road users, the concessionaires that
run tolled roads, and the State agencies that audit them. Penalties for misconduct apply only when the
State audits, and auditing a defecting operator is costly to the State, more so when the State is also
the operator's lender. The model gives two thresholds on the share of users who leave when oversight is
lax: above the first, compliance is stable once it holds; only above the second can it be restored from
widespread misconduct. A rule that rewards compliance while it is rare and punishes defection once it is
common closes the gap between the two. The thresholds are bracketed with 2024 data for the Hanoi–Hai Phong
expressway in Vietnam.

## Contents

| Folder | What it holds |
|---|---|
| `model/` | The model, every computation and check behind the numbers in the paper and the Supporting Information (SI), cached results (`.npy`) and printed results (`.txt`) |
| `figures/` | Scripts that draw the figures of the paper and of the SI at print size; output in `figures/out/` |
| `validation/` | Reproduction of Balabanova et al. (2025), used to validate the pipeline before applying it to the concession game |
| `data/` | The table of seven Vietnamese BOT episodes and the list of all data sources (`SOURCES.md`) |

## Requirements

Python 3.11 or later (tested with Python 3.14) and the packages in `requirements.txt`:

```
pip install -r requirements.txt
```

## Reproducing the results

Run each script from its own folder (`cd model`, then `python <script>`). The cached results are
included, so the figure scripts can be run directly. To recompute everything from scratch, delete the
`.npy` and `.txt` files in `model/` and `figures/` and run the scripts in this order:

| Step | Script | Time |
|---|---|---|
| 1 | `model/verify_symbolic.py` | 1 s |
| 2 | `model/verify_paper_fitness.py` | 2 s |
| 3 | `model/verify_deep.py` | 7 s |
| 4 | `model/verify_block1_block4.py` | 1 s |
| 5 | `model/calibration_2024.py` | 0 s |
| 6 | `model/run_all.py` | 8 min |
| 7 | `model/run_review_checks.py` | 4 min |
| 8 | `model/run_basins.py` | 10 min |
| 9 | `model/run_grid_mp.py` | 43 s |
| 10 | `model/run_mp_checks.py` | 2 min |
| 11 | `model/run_extensions.py` | 9 min |
| 12 | `model/make_figMQ.py` | 16 s |
| 13 | `model/run_calibration_data.py` | 1 s |
| 14 | `model/run_calibration_mc.py` | 1 s |
| 15 | `model/check_pure_reward.py` | 6 min |
| 16 | `figures/fig01_schema.py` | 1 s |
| 17 | `figures/fig02_basins.py` | 5 min |
| 18 | `figures/fig03_05_ode.py` | 19 s |
| 19 | `figures/fig06_cube.py` | 1 s |
| 20 | `figures/fig07_finite.py` | 1 s |
| 21 | `figures/fig08_10_cached.py` | 1 s |
| 22 | `figures/fig13_episodes.py` | 1 s |
| 23 | `figures/si_figs.py` | 2 s |
| 24 | `validation/check_eqs.py` | 1 s |
| 25 | `validation/make_figs.py` | 22 s |

Times are for one run on a laptop (Windows 11, Python 3.14); the full run took about 46 minutes. A fresh run reproduces every cached array and printed result exactly, and every figure of the paper byte for byte.

Figure scripts read the cached results written by the model scripts, so run `model/` first.

## Figures

| Figure | File | Script |
|---|---|---|
| Figure 1 | `figures/out/fig1a_schema.png`, `figures/out/fig1b_replicator.png` | `figures/fig01_schema.py` |
| Figure 2 | `figures/out/fig2_basins.png` | `figures/fig02_basins.py` (data from `model/run_all.py`) |
| Figure 3 | `figures/out/fig8_rule.png` | `figures/fig08_10_cached.py` (data from `model/make_figMQ.py`, `model/run_basins.py`) |
| Figure 4 | `model/fig6_calibration_data.png` | `model/run_calibration_data.py` |
| Figure S1 | `figures/out/fig3_reward.png` | `figures/fig03_05_ode.py` |
| Figure S2 | `figures/out/fig4_lending.png` | `figures/fig03_05_ode.py` |
| Figure S3 | `figures/out/fig5_regimes.png` | `figures/fig03_05_ode.py` |
| Figure S4 | `figures/out/figSI_plane.png` | `figures/si_figs.py` |
| Figure S5 | `figures/out/fig6_cube.png` | `figures/fig06_cube.py` |
| Figure S6 | `figures/out/fig7_finite.png` | `figures/fig07_finite.py` (data from `model/run_all.py`) |
| Figure S7 | `figures/out/figSI_grid.png` | `figures/si_figs.py` (data from `model/run_grid_mp.py`) |
| Figure S8 | `figures/out/figSI_rule.png` | `figures/si_figs.py` (data from `model/run_review_checks.py`) |
| Figure S9 | `figures/out/fig9_users.png` | `figures/fig08_10_cached.py` (data from `model/run_extensions.py`) |
| Figure S10 | `figures/out/fig10_welfare.png` | `figures/fig08_10_cached.py` (data from `model/run_extensions.py`) |
| Figure S11 | `model/figS_calib_mc.png` | `model/run_calibration_mc.py` |
| Figure S12 | `figures/out/fig13_episodes.png` | `figures/fig13_episodes.py` (data in `data/episodes.csv`) |

## Where the numbers come from

| Result in the paper | Script | Where it is printed |
|---|---|---|
| Formulas: fitness functions, fitness differences, thresholds phi*, phi**, phi_s, phi_dagger | `model/verify_symbolic.py`, `model/verify_paper_fitness.py` | console |
| Share of initial states that reach the target (Figure 2) | `model/run_all.py` | `model/figH_basin.npy` |
| Basins under the state dependent rule; phi_dagger = 0.367 at y0 = 1/3 | `model/run_basins.py` | `model/basins.txt` |
| Switching points in finite populations (Z = 50, beta = 0.5): 0.31 punishment, 0.80 pure reward, 0.11 state dependent rule | `model/run_mp_checks.py` | `model/mp_checks.txt` |
| Switching point over population size and selection strength | `model/run_grid_mp.py` | `model/grid_mp.txt` |
| Robustness: gradual rule, reward size, retained fines (kappa), unequal population sizes | `model/run_review_checks.py`, `model/run_mp_checks.py` | `model/review_checks.txt`, `model/mp_checks.txt` |
| Users who observe performance (phi_q^c = 0.79); social welfare of the mechanisms | `model/run_extensions.py` | `model/extensions.txt` |
| Pure reward: cycles for 0.1 < phi < 0.9 | `model/check_pure_reward.py` | console |
| Case study: b_R - b_R5, lambda, phi_emp and the threshold brackets | `model/calibration_2024.py`, `model/run_calibration_data.py`, `model/run_calibration_mc.py` | console |

## Checks

- `model/verify_symbolic.py`: symbolic checks (sympy) of every formula in the paper and the SI.
- `model/verify_paper_fitness.py`, `model/verify_block1_block4.py`: the fitness functions rebuilt from
  the payoff table, compared with the closed forms; the Jacobian at the eight vertices.
- `model/verify_deep.py`: an independent numerical check on random parameter draws, full eigenvalues
  and ODE runs.
- `validation/`: the pipeline reproduces Figures 2 and 3 of Balabanova et al. (2025), and
  `validation/check_eqs.py` confirms their equations (11a–d) and (12a–e) from their payoff tables.

Finite population results use 60 significant digits (mpmath), because at large selection strength the
fixation probabilities fall below the range of double precision.

## Data

- `model/calibration_2024.py`: every number of the case study, each from one official source for 2024.
- `data/episodes.csv`: the seven episodes of toll avoidance or protest (SI, Section 9).
- `data/SOURCES.md`: all sources with links and access dates, and SHA-256 checksums of the source
  documents. The documents themselves (official spreadsheet, audited accounts, tariff notice) are not
  redistributed; download them from the links.

## Licence

The code is released under the MIT License (`LICENSE`). The data in `data/` are compiled from the public
sources listed in `data/SOURCES.md`.

## Citation

If you use this code, please cite the paper:

```
@unpublished{nguyen2026punishment,
  title  = {When Punishment Is Not Enough: Incentive Design for Compliance in Highway
            Public-Private Partnerships},
  author = {Nguyen, Toan Thang and Han, The Anh and Perc, Matja{\v{z}} and Nguyen Ngoc, Doanh},
  year   = {2026},
  note   = {Preprint}
}
```

Contact: The Anh Han (T.Han@tees.ac.uk), Matjaž Perc (matjaz.perc@gmail.com),
Doanh Nguyen Ngoc (doanhnn@vinuni.edu.vn).
