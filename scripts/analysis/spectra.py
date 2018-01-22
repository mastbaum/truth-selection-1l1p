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

def plot_spectra(numu_xvals, inc_numu_spectra, nue_xvals, inc_nue_spectra,
                  nue_nue_spectra, signal_nue_spectra, cov,
                  pot=5e19, eff_m=0.3, eff_e=0.3,
                  truth_eff_m=TRUTH_EFF_M, truth_eff_e=TRUTH_EFF_E, truth_eff_s=TRUTH_EFF_S):
    '''Plot the signal and background spectra.

    :param pot: Exposure (POT)
    :param eff_m: 1m1p efficiency
    :param eff_e: 1e1p efficiency
    '''
    sfb = pot / 5e19
    sfs = pot / 5e19

    scaled_inc_numu = sfb * inc_numu_spectra * eff_m / truth_eff_m
    scaled_inc_nue = sfb * inc_nue_spectra * eff_e / truth_eff_e 
    scaled_nue_nue = sfb * nue_nue_spectra * eff_e / truth_eff_e
    scaled_signal_nue = sfs * signal_nue_spectra * eff_e / truth_eff_s

    scaled_background_nue = scaled_inc_nue + scaled_nue_nue
    scaled_background_numu = scaled_inc_numu

    eb = np.hstack((scaled_inc_numu, scaled_background_nue))
    covr = cov * np.einsum('i,j', eb, eb) #+ np.diag(eb)

    fig, (ax1, ax2) = plt.subplots(1, 2)

    # 1mu1p
    poisi = np.array(scipy.stats.poisson.interval(mu=scaled_background_numu, alpha=0.682))
    stat = np.array([scaled_background_numu - poisi[0], poisi[1] - scaled_background_numu])
    syst = np.sqrt(np.diag(covr[:len(numu_xvals),:len(numu_xvals)]))
    tot  = [np.sqrt(np.square(stat[0]) + np.square(syst)),
            np.sqrt(np.square(stat[1]) + np.square(syst))]

    ax1.bar(numu_xvals, scaled_background_numu, 
        width=np.diff(numu_xvals)[0], color='gray', label='$1\mu1p$')

    ax1.errorbar(numu_xvals, scaled_background_numu, stat, fmt=',',
                 color='black', capsize=3, label='Error (stat)')
    ax1.errorbar(numu_xvals, scaled_background_numu, tot,
                 fmt=',', color='black', label='Error (total)')

    ax1.set_xlabel('$E_\\nu^{CCQE}$ (MeV)')
    ax1.set_ylabel('Events/%1.1f MeV bin/%1.1E POT' % (np.diff(numu_xvals)[0], pot))
    ax1.set_xlim(numu_xvals[0]-np.diff(numu_xvals)[0]/2, numu_xvals[-1]+np.diff(numu_xvals)[0]/2)

    ax1.set_ylim(0, 1.35 * np.max(scaled_background_numu))
    ax1.legend(loc='upper right')
    
    # 1e1p + LEE signal
    poisi = np.array(scipy.stats.poisson.interval(mu=scaled_background_nue, alpha=0.682))
    stat = [scaled_background_nue - poisi[0], poisi[1] - scaled_background_nue]
    syst = np.sqrt(np.diag(covr[len(numu_xvals):,len(numu_xvals):]))
    tot  = [np.sqrt(np.square(stat[0]) + np.square(syst)),
            np.sqrt(np.square(stat[1]) + np.square(syst))]

    ax2.bar(nue_xvals, scaled_background_nue+scaled_signal_nue, 
        width=np.diff(nue_xvals)[0], color='darkred', label='LEE Signal')

    ax2.bar(nue_xvals, scaled_background_nue, 
        width=np.diff(nue_xvals)[0], color='gray', label='Baseline $1e1p$')

    ax2.bar(nue_xvals, scaled_inc_nue,
        width=np.diff(nue_xvals)[0], color='black', label='Mis-ID $1e1p$')

    ax2.errorbar(nue_xvals, scaled_background_nue, stat, fmt=',',
                 color='black', capsize=3, label='Error (stat)')
    ax2.errorbar(nue_xvals, scaled_background_nue, tot,
                 fmt=',', color='black', label='Error (total)')

    import matplotlib.ticker as ticker
    ax2.set_xlabel('$E_\\nu^{CCQE}$ (MeV)')

    ax2.set_xlim(nue_xvals[0]-np.diff(nue_xvals)[0]/2, nue_xvals[-1]+np.diff(nue_xvals)[0]/2)
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(500))
    ax2.set_ylim(0, 1.1 * np.max(scaled_background_nue+scaled_signal_nue))
    ax2.legend(loc='upper right')
    plt.tight_layout()
    plt.show()

def load(background_fname, inclusive_fname, nue_fname, signal_fname):
    '''Load histograms and covariance matrices from disk.

    Nominally these are the output from TSCovariance.

    :param em_file: File with e and mu samples
    :param signal_file: File the the signal sample
    :returns: Tuple of histograms, cov
    '''
    import ROOT

    inclusive_file = ROOT.TFile(inclusive_fname)
    background_file = ROOT.TFile(background_fname)
    nue_file = ROOT.TFile(nue_fname)
    signal_file = ROOT.TFile(signal_fname)

    inclusive_nue = inclusive_file.Get('enu_nue')
    inclusive_numu = inclusive_file.Get('enu_numu')
    nue_nue = nue_file.Get('enu_nue')
    signal = signal_file.Get('hg')
    covariance = background_file.Get('cov')

    # Load raw histograms and covariance matrix

    n_bins_total = inclusive_nue.GetNbinsX() + inclusive_numu.GetNbinsX()
    n_bins_nue = inclusive_nue.GetNbinsX()
    n_bins_numu = inclusive_numu.GetNbinsX()

    cov = np.empty(shape=(n_bins_total, n_bins_total))
    for i in range(n_bins_total):
        for j in range(n_bins_total):
            cov[i][j] = covariance.GetBinContent(i, j)

    # Select desired bins and recompute the covariances
    gmbins = np.arange(2, 20)
    gebins = np.arange(1, n_bins_nue-2)

    good = np.hstack((gmbins, n_bins_numu+gebins))

    inc_numu_spectra = np.array([inclusive_numu.GetBinContent(i) for i in gmbins])
    inc_nue_spectra = np.array([inclusive_nue.GetBinContent(i) for i in gebins])

    nue_nue_spectra = np.array([nue_nue.GetBinContent(i) for i in gebins])

    inc_combined_spectra= np.maximum(1e-40, np.hstack((inc_numu_spectra, inc_nue_spectra)))

    background_combined_spectra = np.maximum(1e-40, \
        np.hstack((inc_numu_spectra, nue_nue_spectra + inc_nue_spectra)))

    cov = cov[good[:,np.newaxis],good] / np.einsum('i,j', background_combined_spectra, background_combined_spectra) 

    numu_xvals = np.array([inclusive_numu.GetBinCenter(int(i)) for i in gmbins])
    nue_xvals = np.array([inclusive_nue.GetBinCenter(int(i)) for i in gebins])

    signal_spectra = np.array([signal.GetBinContent(int(i)) for i in good])
    signal_nue_spectra = signal_spectra[len(numu_xvals):]

    return numu_xvals, inc_numu_spectra, nue_xvals, inc_nue_spectra, nue_nue_spectra, \
        signal_nue_spectra, signal_spectra, cov

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("background_dir")
    parser.add_argument("inclusive_dir")
    parser.add_argument("nue_dir")
    parser.add_argument("signal_dir")
    args = parser.parse_args()

    numu_xvals, inc_numu_spectra, nue_xvals, inc_nue_spectra, nue_nue_spectra, \
    signal_nue_spectra, signal_spectra, cov = \
        load(args.background_dir + '/cov_all.root', 
             args.inclusive_dir + '/cov_all.root',
             args.nue_dir + '/cov_all.root', 
             args.signal_dir + '/cov_all.root')

    plot_spectra(numu_xvals, inc_numu_spectra, nue_xvals, inc_nue_spectra,
        nue_nue_spectra, signal_nue_spectra, cov, pot=6.6e20)

