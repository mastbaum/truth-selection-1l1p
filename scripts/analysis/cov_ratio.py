'''LEE Significance Calculation.

Several functions to compute chi2 significance of a LEE model, and
produce various plots.

The input files and the set of functions that are run are set in the
"if __name__ == '__main__'" code block.

A. Mastbaum <mastbaum@uchicago.edu>, 2017/09/11
'''

import sys
import numpy as np
from numpy.core.umath_tests import inner1d
import scipy.stats

from matplotlib import rc
import os
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
#rc('text', usetex=True)
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['xtick.labelsize'] = 16
matplotlib.rcParams['ytick.labelsize'] = 16
matplotlib.rcParams['axes.titlesize'] = 16
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['legend.fontsize'] = 12

import argparse
import time

TRUTH_EFF_M = 0.379742
TRUTH_EFF_E = 0.782320
TRUTH_EFF_S = 0.970326 

def load(em_file, signal_file):
    '''Load histograms and covariance matrices from disk.

    Nominally these are the output from TSCovariance.

    :param em_file: File with e and mu samples
    :param signal_file: File the the signal sample
    :returns: Tuple of histograms, cov
    '''
    import ROOT

    fi = ROOT.TFile(em_file)
    fs = ROOT.TFile(signal_file)

    henu_e = fi.Get('enu_nue')
    henu_m = fi.Get('enu_numu')
    henu_s = fs.Get('hg')
    hcov = fi.Get('cov')

    # Load raw histograms and covariance matrix
    ebins = range(1, henu_e.GetNbinsX()+1)
    mbins = range(1, henu_m.GetNbinsX()+1)
    em = np.array([henu_m.GetBinContent(i) for i in ebins])
    ee = np.array([henu_e.GetBinContent(i) for i in mbins])

    eb = np.maximum(1e-40, np.hstack((em, ee)))
    cov = np.empty(shape=(len(eb), len(eb)), dtype=eb.dtype)
    for i in range(len(eb)):
        for j in range(len(eb)):
            cov[i][j] = hcov.GetBinContent(i, j)

    # Select desired bins and recompute the covariances
    #gmbins = np.arange(3, 17)
    gmbins = np.arange(2, 20)
    gebins = np.arange(1, len(ee)-2)

    fcov = cov / np.einsum('i,j', eb, eb)
    good = np.hstack((gmbins, len(mbins)+gebins))
    em = em[gmbins]
    ee = ee[gebins]
    eb = eb[good]
    fcov = fcov[good[:,np.newaxis],good]
    cov = fcov * np.einsum('i,j', eb, eb)

    emx = np.array([henu_m.GetBinCenter(int(i)) for i in gmbins])
    eex = np.array([henu_e.GetBinCenter(int(i)) for i in gebins])

    #sscale = 1.0 / 10  # Arbitrary scaling to get the signal ish-right
    sscale = 1.
    es = np.array([henu_s.GetBinContent(int(i)) for i in good]) * sscale
    ess = es[len(emx):]

    return emx, em, eex, ee, es, ess, cov

def plot_cor_ratio(cov_numerator, cov_denomenator, em, outfile='corr.pdf'):
    '''Plot the correlation matrix.

    :param cov: Covariance matrix
    '''
    cor_num = cov_numerator / np.sqrt(np.einsum('ii,jj->ij', cov_numerator, cov_numerator))
    cor_den = cov_denomenator / np.sqrt(np.einsum('ii,jj->ij', cov_denomenator, cov_denomenator))

    cor = np.divide(cor_num, cor_den)
    f, ax = plt.subplots()
    m = ax.matshow(cor, interpolation='nearest', origin='lower', cmap='inferno')
    plt.gca().xaxis.tick_bottom()
    f.colorbar(m)
    ax.set_xlabel('Energy bin index ($\\nu_\\mu,\\nu_e$)')
    ax.set_ylabel('Energy bin index ($\\nu_\\mu,\\nu_e$)')
    ax.axhline(len(em)-0.5, lw=2, color='black')
    ax.axvline(len(em)-0.5, lw=2, color='black')
    plt.tight_layout()
    #if outfile:
    #    plt.savefig(outfile)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("numerator_background_dir")
    parser.add_argument("numerator_signal_dir")
    parser.add_argument("denomenator_background_dir")
    parser.add_argument("denomenator_signal_dir")
    parser.add_argument("--mcmc", action="store_true")
    args = parser.parse_args()

    _, enu_m, _, enu_e, _, _, cov_top = \
        load(args.numerator_background_dir + '/cov_all.root', args.numerator_signal_dir + '/cov_all.root')
    eb = np.hstack((enu_m, enu_e))
    fcov_top = cov_top / np.einsum('i,j', eb, eb)

    _, enu_m, _, enu_e, _, _, cov_bot = \
        load(args.denomenator_background_dir + '/cov_all.root', args.denomenator_signal_dir + '/cov_all.root')
    eb = np.hstack((enu_m, enu_e))
    fcov_bot = cov_bot / np.einsum('i,j', eb, eb)

    # Draw the total correlation matrix
    plot_cor_ratio(fcov_top, fcov_bot, enu_m)

