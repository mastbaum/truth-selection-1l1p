#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TStyle.h"
#include "TSSelection.h"
#include "TSCovariance.h"
#include "TRandom.h"


using namespace std;
namespace galleryfmwk {

std::vector<std::vector<TGraph*> > BinCorrelations(
    TH1D* enu, std::vector<TH1D*> enu_syst) {
  size_t nbins = enu->GetNbinsX();
  size_t nuni = enu_syst.size();

  std::vector<std::vector<TGraph*> > v(nbins);
  for (size_t i=0; i<nbins; i++) {
    v[i] = std::vector<TGraph*>(nbins);
    for (size_t j=0; j<nbins; j++) {
      v[i][j] = new TGraph(nuni);
    }
  }

  for (size_t i=0; i<nbins; i++) {
    for (size_t j=0; j<nbins; j++) {
      for (size_t k=0; k<nuni; k++) {
        v[i][j]->SetPoint(k, enu_syst[k]->GetBinContent(i),
                             enu_syst[k]->GetBinContent(j));
      }
    }
  }

  return v;
}


/******************************************************************************
 ** TSCovariance::EventSample implementation                                 **
 *****************************************************************************/

TSCovariance::EventSample::EventSample(std::string _name,
                                       size_t nbins,
                                       double elo, double ehi,
                                       size_t nweights)
    : name(_name), enu(nullptr), cov(nullptr) {
  enu = new TH1D(("enu_" + name).c_str(),
                 ";CCQE Energy [MeV];Entries per bin",
                 nbins, elo, ehi);
  enu->Sumw2();

  Resize(nweights);
}


TSCovariance::EventSample::~EventSample() {
  delete cov;
  delete enu;
}


TGraphErrors* TSCovariance::EventSample::EnuCollapsed() {
  // Compute the covariance matrix first, if we haven't already
  if (!cov) {
    CovarianceMatrix();
  }

  // x/y values and (symmetric) errors
  const size_t nbins = enu->GetNbinsX();
  double xv[nbins];
  double xe[nbins];
  double yv[nbins];
  double ye[nbins];

  for (size_t i=0; i<nbins; i++) {
    xv[i] = enu->GetBinCenter(i);
    xe[i] = enu->GetBinWidth(i) / 2;
    yv[i] = enu->GetBinContent(i);
    ye[i] = sqrt(cov->GetBinContent(i, i));
  }

  return new TGraphErrors(nbins, xv, yv, xe, ye);
}

void TSCovariance::EventSample::Resize(size_t nweights) {
  enu_syst.clear();
  for (size_t i=0; i<nweights; i++) {
    std::string hname = Form("enu_%s_%zu", name.c_str(), i);
    TH1D* h = (TH1D*) enu->Clone(hname.c_str());
    enu_syst.push_back(h);
  }
}


TH2D* TSCovariance::EventSample::CovarianceMatrix(
    TH1D* nom, std::vector<TH1D*> syst) {
  int nbins = nom->GetNbinsX();

  TH2D* _cov = new TH2D("cov", "", nbins, 0, nbins, nbins, 0, nbins);

  for (int i=1; i<nbins+1; i++) {
    for (int j=1; j<nbins+1; j++) {
      double vij = 0;
      for (size_t k=0; k<syst.size(); k++) {
        double vi = nom->GetBinContent(i) - syst[k]->GetBinContent(i);
        double vj = nom->GetBinContent(j) - syst[k]->GetBinContent(j);
        vij += vi * vj / syst.size();
      }
      _cov->SetBinContent(i, j, vij);
    }
  }

  return _cov;
}


TH2D* TSCovariance::EventSample::CovarianceMatrix() {
  delete cov;
  cov = CovarianceMatrix(enu, enu_syst);
  cov->SetName(("cov_" + name).c_str());
  return cov;
}



TH2D* TSCovariance::EventSample::CorrelationMatrix(TH2D* _cov) {
  TH2D* cor = (TH2D*) _cov->Clone("cor");

  for (int i=1; i<_cov->GetNbinsX()+1; i++) {
    for (int j=1; j<_cov->GetNbinsY()+1; j++) {
      double vij = _cov->GetBinContent(i, j);
      double si = sqrt(_cov->GetBinContent(i, i));
      double sj = sqrt(_cov->GetBinContent(j, j));
      cor->SetBinContent(i, j, vij / (si * sj));
    }
  }

  return cor;
}


TH2D* TSCovariance::EventSample::CorrelationMatrix() {
  // Compute the covariance matrix first, if we haven't already
  if (!cov) {
    CovarianceMatrix();
  }

  TH2D* cor = CorrelationMatrix(cov);
  cor->SetName(("cor_" + name).c_str());
  return cor;
}


/******************************************************************************
 ** TSCovariance implementation                                              **
 *****************************************************************************/

TSCovariance::TSCovariance() : fScaleFactorE(1.0), fScaleFactorMu(1.0),
                               fSeed(0) {}


void TSCovariance::SetEEfficiencyBins(std::vector<float> bins) {
   _energy_efficiency_bins = bins;
}
void TSCovariance::SetAEfficiencyBins(std::vector<float> bins) {
   _angle_efficiency_bins = bins;
}
void TSCovariance::SetEfficiencies(std::vector<float> eff) {
  assert((_energy_efficiency_bins.size() - 1) * (_angle_efficiency_bins.size() - 1) == _efficiencies.size());
  _efficiencies = eff;
}

float TSCovariance::GetEffWeight(float energy, float angle) {
  if (_energy_efficiency_bins.size() == 0) {
    return 1.;
  }

  size_t i;
  for (i = 0; i < _energy_efficiency_bins.size()-1 ; i ++) { 
    if (energy > _energy_efficiency_bins[i] && energy < _energy_efficiency_bins[i]) {
        break;
    }
  } 
  assert(i != _energy_efficiency_bins.size()-1 );
  size_t j;
  for (j = 0; j < _angle_efficiency_bins.size()-1 ; j ++) { 
    if (angle > _angle_efficiency_bins[i] && angle < _angle_efficiency_bins[i]) {
        break;
    }
  } 
  assert(j != _angle_efficiency_bins.size()-1);
  return _efficiencies[j * (_angle_efficiency_bins.size() -1) + i];
}

void TSCovariance::AddWeight(std::string w) {
  use_weights.insert(w);
}


void TSCovariance::init() {
  assert(fInputFile != "" && fOutputFile != "");

  fFile = TFile::Open(fOutputFile.c_str(), "recreate");
  assert(fFile);

  samples.push_back(new EventSample("numu"));
  samples.push_back(new EventSample("nue"));

  // use CCQE energy by default
  _use_ccqe = true;


  std::cout << "TSCovariance: Initialized. Weights: ";
  for (auto it : use_weights) {
    std::cout << it << " ";
  }
  std::cout << std::endl;

  std::cout << "TSCovariance: Writing output to " << fOutputFile << std::endl;
  _efficiencies = std::vector<float>();
  _energy_efficiency_bins = std::vector<float>();
  _angle_efficiency_bins = std::vector<float>();

  gRandom->SetSeed(fSeed);
}

void TSCovariance::setUseCCQE(bool b) {
  _use_ccqe = b;
} 

void TSCovariance::analyze() {
  // Grab the MCTruth information
  TFile f(fInputFile.c_str());
  TTree* _tree = (TTree*) f.Get("data");
  assert(_tree && _tree->GetEntries() > 0);

  galleryfmwk::TSSelection::OutputData _data;
  _tree->SetBranchAddress("eccqe", &_data.eccqe);
  _tree->SetBranchAddress("lpid", &_data.lpid);
  _tree->SetBranchAddress("lpdg", &_data.lpdg);
  _tree->SetBranchAddress("nupdg", &_data.nupdg);
  _tree->SetBranchAddress("ccnc", &_data.ccnc);
  _tree->SetBranchAddress("bnbweight", &_data.bnbweight);
  _tree->SetBranchAddress("weights", &_data.weights);
  //_tree->SetBranchAddress("eps", &_data.eps);
  _tree->SetBranchAddress("elep", &_data.elep);
  _tree->SetBranchAddress("thetalep", &_data.thetalep);

  // Event loop
  for (long k=0; k<_tree->GetEntries(); k++) {
    _tree->GetEntry(k);

    // Iterate through all the weighting functions to compute a set of
    // total weights for this event. mcWeight is a mapping from reweighting
    // function name to a vector of weights for each "universe."
    std::vector<double> weights;
    size_t wmin = 1000000;
    for (auto const& it : *_data.weights) {
      if (use_weights.find(it.first) == use_weights.end()) { continue; }
      if (it.second.size() < wmin) {
        wmin = it.second.size();
      }

      assert(wmin < 1000000);
      weights.resize(wmin, 1.0);

      // Compute universe-wise product of all requsted weights
      if (use_weights.find("*") != use_weights.end() ||
          use_weights.find(it.first) != use_weights.end()) {
        for (size_t i=0; i<weights.size(); i++) {
          weights[i] *= it.second[i];
        }
      }
    }

    // The observable
    double nuEnergy = _use_ccqe ? _data.eccqe : 
       std::accumulate(_data.eps.begin(), _data.eps.end(), 0.0) + _data.elep;

    // Determine which event sample this event corresponds to, based on
    // lepton PID
    EventSample* sample = nullptr;
    double fs = 1.0;  // Scale factor

    // Scale for CCnue or inclusive sample
    if (abs(_data.nupdg) == 12 && _data.ccnc == 0) {
      fs = fScaleFactorE;
    }
    else {
      fs = fScaleFactorMu;
    }

    if (_data.lpid == 13) {
      sample = samples[0];
    }
    else if (_data.lpid == 11) {
      sample = samples[1];
    }
    else {
      std::cout << "Unknown lepton PID " << _data.lpid << std::endl;
      continue;
    }

    fs *= _data.bnbweight;  // Apply BNB correction weight
    // apply efficiency weight (binned for now by energy and angle)
    fs *= GetEffWeight(nuEnergy, _data.thetalep); 

    // Fill histograms for this event sample
    if (sample->enu_syst.empty()) {
      sample->Resize(weights.size());
    }
    else {
      assert(sample->enu_syst.size() == weights.size());
    }

    // CV histogram
    sample->enu->Fill(nuEnergy, fs);

    // Systematics histograms with weights
    for (size_t i=0; i<weights.size(); i++) {
      sample->enu_syst[i]->Fill(nuEnergy, weights[i] * fs);
    }
  }

  /////////////////////////////////////////////////////////
  // Output
  size_t total_bins = 0;

  fFile->cd();

  // Write out sample=-wise distributions
  for (size_t i=0; i<samples.size(); i++) {
    samples[i]->enu->Write();
    total_bins += samples[i]->enu->GetNbinsX();

    TH2D* cov = samples[i]->CovarianceMatrix();
    cov->Write();

    TH2D* cor = samples[i]->CorrelationMatrix();
    cor->Write();

    TGraphErrors* g = samples[i]->EnuCollapsed();
    g->SetName(("err_" + samples[i]->name).c_str());
    g->Write();
  }

  // Global (sample-to-sample) distributions
  // Build glued-together energy spectra for the nominal and each systematics
  // universe, and feed into the correlation matrix calculator.
  total_bins -= samples.size();
  TH1D hg("hg", ";E_{#nu};Entries per bin", total_bins, 0, total_bins);
  hg.Sumw2();
  std::vector<TH1D*> hgsys;
  for (size_t i=0; i<samples[0]->enu_syst.size(); i++) {
    hgsys.push_back(new TH1D(Form("hg%zu", i), "", total_bins, 0, total_bins));
    hgsys[i]->Sumw2();
  }
  size_t ibin = 0;
  for (size_t i=0; i<samples.size(); i++) {
    for(int j=1; j<samples[i]->enu->GetNbinsX()+1; j++) {
      hg.SetBinContent(ibin, samples[i]->enu->GetBinContent(j));
      hg.SetBinError(ibin, samples[i]->enu->GetBinError(j));
      if (hgsys.size() != samples[i]->enu_syst.size()) {
        ibin++;
        continue;
      }
      for (size_t k=0; k<hgsys.size(); k++) {
        hgsys[k]->SetBinContent(ibin, samples[i]->enu_syst[k]->GetBinContent(j));
      }
      ibin++;
    }
  }

  hg.Write();

  TH2D* gcov = EventSample::CovarianceMatrix(&hg, hgsys);
  gcov->Write();

  TH2D* gcor = EventSample::CorrelationMatrix(gcov);
  gcor->Write();

  // Write plots of bin correlations for each pair of bins
  std::vector<std::vector<TGraph*> > vv = BinCorrelations(&hg, hgsys);
  for (long i=0; i<hg.GetNbinsX(); i++) {
    for (long j=0; j<hg.GetNbinsX(); j++) {
      char nname[50];
      snprintf(nname, 50, "gc_%02lu_%02lu", i, j);
      vv[i][j]->SetName(nname);
      vv[i][j]->Write();
    }
  }
}

}  // namespace galleryfmwk

