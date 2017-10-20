'''Produce plots of covariances.

Makes PDF format plots from the output of `cov.py` (files which contain
covariance matrices and spectra with error bands).

Usage:

  $ python plotcov.py

This script should be modified to use the desired combinations of
systematic parameters, and the input filenames.

A. Mastbaum <mastbaum@uchicago.edu>, 2017/09
'''

import ROOT
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# Names of the systematics to plot
wnames = [
    #'expskin_FluxUnisim',
    #'genie_AGKYpT_Genie',
    #'genie_AGKYxF_Genie',
    #'genie_DISAth_Genie',
    #'genie_DISBth_Genie',
    #'genie_DISCv1u_Genie',
    #'genie_DISCv2u_Genie',
    #'genie_FermiGasModelKf_Genie',
    #'genie_FermiGasModelSf_Genie',
    #'genie_FormZone_Genie',
    #'genie_IntraNukeNabs_Genie',
    #'genie_IntraNukeNcex_Genie',
    #'genie_IntraNukeNel_Genie',
    #'genie_IntraNukeNinel_Genie',
    #'genie_IntraNukeNmfp_Genie',
    #'genie_IntraNukeNpi_Genie',
    #'genie_IntraNukePIabs_Genie',
    #'genie_IntraNukePIcex_Genie',
    #'genie_IntraNukePIel_Genie',
    #'genie_IntraNukePIinel_Genie',
    #'genie_IntraNukePImfp_Genie',
    #'genie_IntraNukePIpi_Genie',
    #'genie_NC_Genie',
    #'genie_NonResRvbarp1pi_Genie',
    #'genie_NonResRvbarp2pi_Genie',
    #'genie_NonResRvp1pi_Genie',
    #'genie_NonResRvp2pi_Genie',
    #'genie_ResDecayEta_Genie',
    #'genie_ResDecayGamma_Genie',
    #'genie_ResDecayTheta_Genie',
    #'genie_ccresAxial_Genie',
    #'genie_ccresVector_Genie',
    #'genie_cohMA_Genie',
    #'genie_cohR0_Genie',
    #'genie_ncelAxial_Genie',
    #'genie_ncelEta_Genie',
    #'genie_ncresAxial_Genie',
    #'genie_ncresVector_Genie',
    #'genie_qema_Genie',
    #'genie_qevec_Genie',
    #'horncurrent_FluxUnisim',
    #'kminus_PrimaryHadronNormalization',
    #'kplus_PrimaryHadronFeynmanScaling',
    #'kzero_PrimaryHadronSanfordWang',
    #'nucleoninexsec_FluxUnisim',
    #'nucleonqexsec_FluxUnisim',
    #'nucleontotxsec_FluxUnisim',
    #'piminus_PrimaryHadronSWCentralSplineVariation',
    #'pioninexsec_FluxUnisim',
    #'pionqexsec_FluxUnisim',
    #'piontotxsec_FluxUnisim',
    #'piplus_PrimaryHadronSWCentralSplineVariation',
    #'xsr_scc_Fv3_XSecRatio',
    #'xsr_scc_Fa3_XSecRatio',
    'genie',
    'flux',
    'xsr',
    'all',
]

for w in wnames:
    print w
    f = ROOT.TFile('cov_em/cov_%s.root' % w)
    c1 = ROOT.TCanvas()

    he = f.Get('err_nue')
    he.SetFillColor(ROOT.kRed-7)
    he.SetTitle(';CCQE Energy (MeV);Events')
    he.Draw('a2')
    c1.SaveAs('figures/nue_%s.pdf' % w)

    hm = f.Get('err_numu')
    hm.SetFillColor(ROOT.kRed-7)
    hm.SetTitle(';CCQE Energy (MeV);Events')
    hm.Draw('a2')
    c1.SaveAs('figures/num_%s.pdf' % w)

    cor = f.Get('cor')
    cor.Draw('colz')
    c1.SetRightMargin(0.18)
    c1.SaveAs('figures/cor_%s.pdf' % w)

    cov = f.Get('cov')
    fcov = cov.Clone('fcov')
    fcov.Reset()
    hg = f.Get('hg')

    g = []
    diag = []
    for i in range(0, hg.GetNbinsX()):
        g.append(hg.GetBinContent(i))
        for j in range(0, hg.GetNbinsX()):
            if i == j: diag.append(cov.GetBinContent(i, j))
            d = hg.GetBinContent(i) * hg.GetBinContent(j)
            print('%i %i %f %f' % (i, j, d, cov.GetBinContent(i,j)))
            if d > 1e-5:
                fcov.SetBinContent(i, j, cov.GetBinContent(i, j) / d)

    fcov.Draw('colz')
    c1.SetRightMargin(0.18)
    c1.SaveAs('figures/fcov_%s.pdf' % w)

    g = np.array(g)
    diag = np.array(diag)
    print np.sqrt(diag) / g

