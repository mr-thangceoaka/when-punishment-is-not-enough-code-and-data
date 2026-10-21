"""Common style for the print-size figures.

Figures are drawn at the size at which they are printed (text width of the paper or of the SI), so
they are inserted at scale 1 and all text is at least 8 pt. Figures that rely on cached results
read the .npy files written by the scripts in ../model; quick computations are rerun.
Output: figures/out/*.png at 300 dpi.

Code for "When Punishment Is Not Enough: Incentive Design for Compliance in Highway
Public-Private Partnerships" (T. T. Nguyen, T. A. Han, M. Perc, D. Nguyen Ngoc). See README.md.
"""
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RERUN = os.path.normpath(os.path.join(HERE, '..', 'model'))
DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
OUT = os.path.join(HERE, 'out')
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, RERUN)

TW = 360.0 / 72.27
TW_SI = 472.03 / 72.27

FS = 8.5


def style():
    plt.rcParams.update({
        'font.size': FS, 'axes.labelsize': 9, 'axes.titlesize': 9,
        'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
        'mathtext.fontset': 'dejavusans',
        'lines.linewidth': 1.4, 'axes.linewidth': 0.7,
        'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
        'xtick.major.size': 3, 'ytick.major.size': 3,
        'xtick.major.pad': 2, 'ytick.major.pad': 2, 'axes.labelpad': 2,
        'legend.framealpha': 0.95, 'legend.borderpad': 0.3, 'legend.handlelength': 1.6,
        'legend.labelspacing': 0.25, 'legend.handletextpad': 0.4,
        'figure.dpi': 100, 'savefig.dpi': 300,
    })


def rr(name):
    """Path of a result file written by the scripts in ../model."""
    return os.path.join(RERUN, name)


def save(fig, name):
    """Save at the exact figure size (no cropping), so the figure is inserted at scale 1."""
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=300, facecolor='white')
    plt.close(fig)
    w, h = fig.get_size_inches()
    print(f'saved {name}: {w:.2f} x {h:.2f} in')
