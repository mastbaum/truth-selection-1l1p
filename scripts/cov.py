import sys
from glob import glob
from ROOT import galleryfmwk

#wnames = [
#    'expskin_FluxUnisim',
#    'genie_AGKYpT_Genie',
#    'genie_AGKYxF_Genie',
#    'genie_DISAth_Genie',
#    'genie_DISBth_Genie',
#    'genie_DISCv1u_Genie',
#    'genie_DISCv2u_Genie',
#    'genie_FermiGasModelKf_Genie',
#    'genie_FermiGasModelSf_Genie',
#    'genie_FormZone_Genie',
#    'genie_IntraNukeNabs_Genie',
#    'genie_IntraNukeNcex_Genie',
#    'genie_IntraNukeNel_Genie',
#    'genie_IntraNukeNinel_Genie',
#    'genie_IntraNukeNmfp_Genie',
#    'genie_IntraNukeNpi_Genie',
#    'genie_IntraNukePIabs_Genie',
#    'genie_IntraNukePIcex_Genie',
#    'genie_IntraNukePIel_Genie',
#    'genie_IntraNukePIinel_Genie',
#    'genie_IntraNukePImfp_Genie',
#    'genie_IntraNukePIpi_Genie',
#    'genie_NC_Genie',
#    'genie_NonResRvbarp1pi_Genie',
#    'genie_NonResRvbarp2pi_Genie',
#    'genie_NonResRvp1pi_Genie',
#    'genie_NonResRvp2pi_Genie',
#    'genie_ResDecayEta_Genie',
#    'genie_ResDecayGamma_Genie',
#    'genie_ResDecayTheta_Genie',
#    'genie_ccresAxial_Genie',
#    'genie_ccresVector_Genie',
#    'genie_cohMA_Genie',
#    'genie_cohR0_Genie',
#    'genie_ncelAxial_Genie',
#    'genie_ncelEta_Genie',
#    'genie_ncresAxial_Genie',
#    'genie_ncresVector_Genie',
#    'genie_qema_Genie',
#    'genie_qevec_Genie',
#    'horncurrent_FluxUnisim',
#    'kminus_PrimaryHadronNormalization',
#    'kplus_PrimaryHadronFeynmanScaling',
#    'kzero_PrimaryHadronSanfordWang',
#    'nucleoninexsec_FluxUnisim',
#    'nucleonqexsec_FluxUnisim',
#    'nucleontotxsec_FluxUnisim',
#    'piminus_PrimaryHadronSWCentralSplineVariation',
#    'pioninexsec_FluxUnisim',
#    'pionqexsec_FluxUnisim',
#    'piontotxsec_FluxUnisim',
#    'piplus_PrimaryHadronSWCentralSplineVariation',
#    'xsr_scc_Fv3_XSecRatio',
#    'xsr_scc_Fa3_XSecRatio',
#]

wg = {
    'genie': [
        'genie_AGKYpT_Genie',
        'genie_AGKYxF_Genie',
        'genie_DISAth_Genie',
        'genie_DISBth_Genie',
        'genie_DISCv1u_Genie',
        'genie_DISCv2u_Genie',
        'genie_FermiGasModelKf_Genie',
        'genie_FermiGasModelSf_Genie',
        'genie_FormZone_Genie',
        'genie_IntraNukeNabs_Genie',
        'genie_IntraNukeNcex_Genie',
        'genie_IntraNukeNel_Genie',
        'genie_IntraNukeNinel_Genie',
        'genie_IntraNukeNmfp_Genie',
        'genie_IntraNukeNpi_Genie',
        'genie_IntraNukePIabs_Genie',
        'genie_IntraNukePIcex_Genie',
        'genie_IntraNukePIel_Genie',
        'genie_IntraNukePIinel_Genie',
        'genie_IntraNukePImfp_Genie',
        'genie_IntraNukePIpi_Genie',
        'genie_NC_Genie',
        'genie_NonResRvbarp1pi_Genie',
        'genie_NonResRvbarp2pi_Genie',
        'genie_NonResRvp1pi_Genie',
        'genie_NonResRvp2pi_Genie',
        'genie_ResDecayEta_Genie',
        'genie_ResDecayGamma_Genie',
        'genie_ResDecayTheta_Genie',
        'genie_ccresAxial_Genie',
        'genie_ccresVector_Genie',
        'genie_cohMA_Genie',
        'genie_cohR0_Genie',
        'genie_ncelAxial_Genie',
        'genie_ncelEta_Genie',
        'genie_ncresAxial_Genie',
        'genie_ncresVector_Genie',
        'genie_qema_Genie',
        'genie_qevec_Genie',
    ],
    'flux': [
        'expskin_FluxUnisim',
        'horncurrent_FluxUnisim',
        'kminus_PrimaryHadronNormalization',
        'kplus_PrimaryHadronFeynmanScaling',
        'kzero_PrimaryHadronSanfordWang',
        'nucleoninexsec_FluxUnisim',
        'nucleonqexsec_FluxUnisim',
        'nucleontotxsec_FluxUnisim',
        'piminus_PrimaryHadronSWCentralSplineVariation',
        'pioninexsec_FluxUnisim',
        'pionqexsec_FluxUnisim',
        'piontotxsec_FluxUnisim',
        'piplus_PrimaryHadronSWCentralSplineVariation',
    ],
    'xsr': [
        'xsr_scc_Fv3_XSecRatio',
        'xsr_scc_Fa3_XSecRatio',
    ],
}

wg['all'] = []
for v in wg.values():
    wg['all'] += v

sfe = 5.0e19 / 8.47443e+21
sfm = 5.0e19 / 2.81456e+19

sfs = 5.0e19 / 1.86828e+21
sfe=sfs

#m = galleryfmwk.TSCovariance()
#m.SetInputFile(sys.argv[1])
#m.SetScaleFactorE(sfe)
#m.SetScaleFactorMu(sfm)
#for w in wnames:
#    m.AddWeight(w)
#m.SetOutputFile('cov_all.root')
#m.init()
#m.analyze()
#del m

for w, ws in wg.items():
    print w
    m = galleryfmwk.TSCovariance()
    m.SetInputFile(sys.argv[1])
    m.SetScaleFactorE(sfe)
    m.SetScaleFactorMu(sfm)
    for ww in ws:
        m.AddWeight(ww)
    m.SetOutputFile('cov_%s.root' % w)
    m.init()
    m.analyze()
    del m

for w in wg['all']: #wnames:
    print w
    m = galleryfmwk.TSCovariance()
    m.SetInputFile(sys.argv[1])
    m.SetScaleFactorE(sfe)
    m.SetScaleFactorMu(sfm)
    m.AddWeight(w)
    m.SetOutputFile('cov_%s.root' % w)
    m.init()
    m.analyze()
    del m

