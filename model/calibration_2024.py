"""Calibration constants for the case study (Hanoi-Hai Phong expressway, base year 2024).

Each quantity comes from one official source (links and access dates in data/SOURCES.md and in the
bibliography of the SI):
  [DRVN]     Vietnam Road Administration, toll revenue and traffic of the BOT road projects,
             January to December 2024: rows of the expressway station and of the two QL5 stations.
             Revenue includes 8% VAT. The table lists 53 projects and 58 stations.
  [FS]       VIDIFI, audited consolidated financial statements 2024: interest expense (income
             statement), loans from VDB (note V.13), VDB share of equity (note V.19), toll revenue
             excluding VAT (note VI.1).
  [NOTICE]   VIDIFI tariff notice effective 29 December 2023: full-trip expressway fares and
             single-trip fares per QL5 station, five vehicle classes.
  [BDT2016]  Bao Dau thau, 17 October 2016: operator's count of about 30,000 vehicles per day
             through the QL5 stations and 40,000 to 45,000 on the road.
  [TP2019], [NDT2017], [TCDN2019]: national ranges of toll collection costs (6% to 12% of revenue
             under manual collection; 2.15% to 4.5% under electronic collection).
Cross-check: DRVN revenue = FS revenue x 1.08 for the expressway and for QL5.
Money in million VND unless stated. Run this file to print the derived quantities.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
CT04_REV_2024 = 2_857_454.517
CT04_TRAFFIC_2024 = 21_902_789
QL5_REV_2024 = 994_510.183
QL5_PASSAGES_2024 = 14_168_750
DRVN_N_PROJECTS = 53
DRVN_N_STATIONS = 58
FS_CT04_REV_EXVAT = 2_645_791.220645
FS_QL5_REV_EXVAT = 920_842.762077
FS_INTEREST_2024 = 1_695_587.138699
FS_VDB_DEBT_END2024 = 22_825_889.324512
FS_VDB_SHARE = 0.9697
CT04_FULL_TRIP = (216000, 309000, 411000, 658000, 848000)
QL5_PER_STATION = (44000, 64000, 85000, 137000, 196000)
QL5_THROUGH_2016 = 30_000
QL5_ONROUTE_2016 = (40_000, 45_000)
CR_PRE = (0.06, 0.12)
CR_POST = (0.02, 0.045)
V_ALLOW = (0.0, 0.1)

B_R_VND = CT04_REV_2024 * 1e6 / CT04_TRAFFIC_2024
R_REALISED = (QL5_REV_2024 / QL5_PASSAGES_2024) / (CT04_REV_2024 / CT04_TRAFFIC_2024)
R_CLASS = tuple(2 * q / t for q, t in zip(QL5_PER_STATION, CT04_FULL_TRIP))
R_TARIFF = (min(R_CLASS), max(R_CLASS))
DELTA = (1 - R_TARIFF[1], 1 - R_TARIFF[0])
LAM = FS_INTEREST_2024 / CT04_REV_2024
VLAM = (LAM + V_ALLOW[0], LAM + V_ALLOW[1])
K_OVER_R = FS_VDB_DEBT_END2024 / CT04_REV_2024
PHI_EMP = (1 - QL5_THROUGH_2016 / QL5_ONROUTE_2016[0], 1 - QL5_THROUGH_2016 / QL5_ONROUTE_2016[1])

assert abs(FS_CT04_REV_EXVAT * 1.08 - CT04_REV_2024) < 0.01
assert abs(FS_QL5_REV_EXVAT * 1.08 - QL5_REV_2024) < 0.01

if __name__ == '__main__':
    import sys
    sys.stdout.reconfigure(encoding='utf-8')
    print(f'b_R = {B_R_VND:,.0f} VND/vehicle; realised r = {R_REALISED:.4f}')
    print('r by vehicle class:', ' '.join(f'{r:.4f}' for r in R_CLASS), '-> b_R-b_R5 =', [round(d, 4) for d in DELTA])
    print(f'lambda = {LAM:.4f}; v+lambda = [{VLAM[0]:.4f}, {VLAM[1]:.4f}]; K/R = {K_OVER_R:.3f} years')
    print(f'phi_emp = [{PHI_EMP[0]:.4f}, {PHI_EMP[1]:.4f}]')
