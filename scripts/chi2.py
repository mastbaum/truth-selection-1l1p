import sys
import numpy as np
import ROOT

laura= [
    ROOT.kRed,
    ROOT.TColor.GetColor(153, 153, 153),
    ROOT.TColor.GetColor(75,  110, 188),
    ROOT.TColor.GetColor(110, 157, 100),
    ROOT.TColor.GetColor(121, 93,  136),
    ROOT.TColor.GetColor(255, 153, 4),
    ROOT.TColor.GetColor(153, 102, 51),
    ROOT.TColor.GetColor(0,   204, 204)
]

fi = ROOT.TFile('cov_em/cov_all.root')
fs = ROOT.TFile('cov_sig/cov_all.root')

henu_e = fi.Get('enu_nue')
henu_m = fi.Get('enu_numu')
henu_s = fs.Get('hg')
hcov = fi.Get('cov')

henu_b = henu_e.Clone('henu_b')
henu_b.Add(henu_m)
henu_b.SetFillColor(laura[0])
henu_b.Draw()
henu_t = henu_s.Clone('enu_t')
henu_t.Add(henu_b)
henu_t.SetLineColor(ROOT.kBlack)
henu_t.Draw('same')

brange = range(3,19) + range(24,hcov.GetNbinsX()+1)

enu_m = np.array([henu_m.GetBinContent(i) for i in range(3, 19)])
enu_e = np.array([henu_e.GetBinContent(i) for i in range(1, henu_e.GetNbinsX()+1)])
enu_b = np.hstack((enu_m, enu_e))
nbins = len(enu_b)

print len(brange), len(enu_m), len(enu_e), len(enu_b)

enu_s = np.array([henu_s.GetBinContent(i) for i in brange])

cov = np.empty(shape=(nbins, nbins), dtype=enu_b.dtype)
for i, hi in enumerate(brange):
    for j, hj in enumerate(brange):
        cov[i][j] = hcov.GetBinContent(hi, hj)

print cov

cov += np.diag(enu_b)

icov = np.linalg.inv(cov)
#icov = 1.0/cov

print enu_b
print enu_s

print cov.shape, enu_b.shape, enu_s.shape

chi2 = np.einsum('i,ij,j', enu_b - enu_s, icov, enu_b - enu_s)

print chi2


raw_input()

