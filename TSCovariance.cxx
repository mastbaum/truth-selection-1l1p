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

namespace galleryfmwk {

/******************************************************************************
 ** TSCovariance::EventSample implementation                             **
 *****************************************************************************/

TSCovariance::EventSample::EventSample(std::string _name,
                                            size_t nbins,
                                            double elo, double ehi,
                                            size_t nweights)
    : name(_name), enu(nullptr), cov(nullptr) {
  enu = new TH1D(("enu_" + name).c_str(),
                 ";CCQE Energy [MeV];Entries per bin",
                 nbins, elo, ehi);

  Resize(nweights);
}


TSCovariance::EventSample::~EventSample() {
  delete cov;
  delete enu;
}


TGraphErrors* TSCovariance::EventSample::EnuCollapsed() {
  size_t nbins = enu->GetNbinsX();

  // Compute the mean and standard deviation across universes using
  // Welford's method, cf. Art of Computer Programming (D. Knuth)
  double x[nbins];
  double y[nbins];
  double y0[nbins];
  double s[nbins];
  double s0[nbins];
  double xe[nbins];

  if (enu_syst.empty()) {
    std::cout << "TSCovariance::EventSample::EnuCollapsed: "
              << "No multisims for sample " << name << std::endl;
    return new TGraphErrors();
  }

  for (size_t i=0; i<nbins; i++) {
    x[i] = enu->GetBinCenter(i);
    y[i] = y0[i] = enu_syst[0]->GetBinContent(i+1);
    s[i] = s0[i] = 0.0;
    xe[i] = enu->GetBinWidth(i) / 2;

    for (size_t j=1; j<enu_syst.size(); j++) {
      double v = enu_syst[j]->GetBinContent(i);
      y[i] = y0[i] + (v - y0[i]) / j;
      s[i] = s0[i] + (v - y0[i]) * (v - y[i]);

      y0[i] = y[i];
      s0[i] = s[i];
    }

    s[i] = sqrt(s[i] / enu_syst.size());
  }

  return new TGraphErrors(nbins, x, y, xe, s);
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
      for (size_t k = 0; k<syst.size(); k++) {
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
 ** TSCovariance implementation                                          **
 *****************************************************************************/

TSCovariance::TSCovariance() : fScaleFactorE(1.0), fScaleFactorMu(1.0) {}


TSCovariance::~TSCovariance() {
  ofs.close();
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

  //std::vector<std::string> wv = {
  //  // "*"
  //  //,"bnbcorrection_FluxHist"
  //   "expskin_FluxUnisim"
  //  ,"genie_AGKYpT_Genie"
  //  ,"genie_AGKYxF_Genie"
  //  ,"genie_DISAth_Genie"
  //  ,"genie_DISBth_Genie"
  //  ,"genie_DISCv1u_Genie"
  //  ,"genie_DISCv2u_Genie"
  //  ,"genie_FermiGasModelKf_Genie"
  //  ,"genie_FermiGasModelSf_Genie"
  //  ,"genie_FormZone_Genie"
  //  ,"genie_IntraNukeNabs_Genie"
  //  ,"genie_IntraNukeNcex_Genie"
  //  ,"genie_IntraNukeNel_Genie"
  //  ,"genie_IntraNukeNinel_Genie"
  //  ,"genie_IntraNukeNmfp_Genie"
  //  ,"genie_IntraNukeNpi_Genie"
  //  ,"genie_IntraNukePIabs_Genie"
  //  ,"genie_IntraNukePIcex_Genie"
  //  ,"genie_IntraNukePIel_Genie"
  //  ,"genie_IntraNukePIinel_Genie"
  //  ,"genie_IntraNukePImfp_Genie"
  //  ,"genie_IntraNukePIpi_Genie"
  //  ,"genie_NC_Genie"
  //  ,"genie_NonResRvbarp1pi_Genie"
  //  ,"genie_NonResRvbarp2pi_Genie"
  //  ,"genie_NonResRvp1pi_Genie"
  //  ,"genie_NonResRvp2pi_Genie"
  //  ,"genie_ResDecayEta_Genie"
  //  ,"genie_ResDecayGamma_Genie"
  //  ,"genie_ResDecayTheta_Genie"
  //  ,"genie_ccresAxial_Genie"
  //  ,"genie_ccresVector_Genie"
  //  ,"genie_cohMA_Genie"
  //  ,"genie_cohR0_Genie"
  //  ,"genie_ncelAxial_Genie"
  //  ,"genie_ncelEta_Genie"
  //  ,"genie_ncresAxial_Genie"
  //  ,"genie_ncresVector_Genie"
  //  ,"genie_qema_Genie"
  //  ,"genie_qevec_Genie"
  //  ,"horncurrent_FluxUnisim"
  //  ,"kminus_PrimaryHadronNormalization"
  //  ,"kplus_PrimaryHadronFeynmanScaling"
  //  ,"kzero_PrimaryHadronSanfordWang"
  //  ,"nucleoninexsec_FluxUnisim"
  //  ,"nucleonqexsec_FluxUnisim"
  //  ,"nucleontotxsec_FluxUnisim"
  //  ,"piminus_PrimaryHadronSWCentralSplineVariation"
  //  ,"pioninexsec_FluxUnisim"
  //  ,"pionqexsec_FluxUnisim"
  //  ,"piontotxsec_FluxUnisim"
  //  ,"piplus_PrimaryHadronSWCentralSplineVariation"
  //};

  //// "Fv3technote_XSecRatio", "*"

  //use_weights = std::set<std::string>(wv.begin(), wv.end());

  //std::string out_file = "cov.txt";
  //ofs.open(out_file, std::ofstream::out);

  std::cout << "TSCovariance: Initialized. Weights: ";
  for (auto it : use_weights) {
    std::cout << it << " ";
  }
  std::cout << std::endl;

  std::cout << "TSCovariance: Writing output to " << fOutputFile << std::endl;
}
  

void TSCovariance::analyze() {
  // Grab the MCTruth information
  TFile f(fInputFile.c_str());
  TTree* _tree = (TTree*) f.Get("data");
  assert(_tree && _tree->GetEntries() > 0);

  galleryfmwk::TSSelection::OutputData _data;
  _tree->SetBranchAddress("enu", &_data.enu);
  _tree->SetBranchAddress("q2", &_data.q2);
  _tree->SetBranchAddress("w", &_data.w);
  _tree->SetBranchAddress("int", &_data.int_type);
  _tree->SetBranchAddress("mode", &_data.int_mode);
  _tree->SetBranchAddress("ccnc", &_data.ccnc);
  _tree->SetBranchAddress("eccqe", &_data.eccqe);
  _tree->SetBranchAddress("ep", &_data.ep);
  _tree->SetBranchAddress("ppdg", &_data.ppdg);
  _tree->SetBranchAddress("elep", &_data.elep);
  _tree->SetBranchAddress("lpdg", &_data.lpdg);
  _tree->SetBranchAddress("lpid", &_data.lpid);
  _tree->SetBranchAddress("bnbweight", &_data.bnbweight);
  _tree->SetBranchAddress("dataset", &_data.dataset);
  _tree->SetBranchAddress("weights", &_data.weights);

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

      //std::cout << "TSCovariance: Found weight: " << it.first << std::endl;

      // Compute universe-wise product of all requsted weights
      if (use_weights.find("*") != use_weights.end() ||
          use_weights.find(it.first) != use_weights.end()) {
        for (size_t i=0; i<weights.size(); i++) {
          weights[i] *= it.second[i];
        }
      }
    }

    // Neutrino interaction truth
    double nuEnergy = _data.eccqe;

    // Determine which event sample this event corresponds to
    EventSample* sample = nullptr;
    double fs = 1.0;

    if (_data.lpid == 11) {
      sample = samples[1];
      fs = fScaleFactorE;
    }
    else if (_data.lpid == 13) {
      sample = samples[0];
      fs = fScaleFactorMu;
    }
    else {
      std::cout << "Unknown lepton PDG " << _data.lpid << std::endl;
      continue;
    }
    //std::cout << fs << " ";
    fs *= _data.bnbweight;
    //std::cout << _data.bnbweight << " " << fs << std::endl;

    // Fill histograms for this event sample
    if (sample->enu_syst.empty()) {
      sample->Resize(weights.size());
    }
    else {
      assert(sample->enu_syst.size() == weights.size());
    }

    // Neutrino Energy
    sample->enu->Fill(nuEnergy, fs);

    // Fill weights
    //std::cout << "TSCovariance: Weights: ";
    ofs << _data.lpdg << "\t" << nuEnergy << "\t";
    for (size_t i=0; i<weights.size(); i++) {
      sample->enu_syst[i]->Fill(nuEnergy, weights[i] * fs);
      //std::cout << weights[i] << " ";
      ofs << weights[i] << "\t";
    }
    //std::cout << std::endl;
    ofs << std::endl;
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
  std::vector<TH1D*> hgsys;
  for (size_t i=0; i<samples[0]->enu_syst.size(); i++) {
    hgsys.push_back(new TH1D(Form("hg%zu", i), "", total_bins, 0, total_bins));
  }
  size_t ibin = 0;
  for (size_t i=0; i<samples.size(); i++) {
    for(int j=1; j<samples[i]->enu->GetNbinsX()+1; j++) {
      hg.SetBinContent(ibin, samples[i]->enu->GetBinContent(j));
      if (hgsys.size() != samples[i]->enu_syst.size()) {
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
}

TGraphErrors* EnuCollapsed2(std::vector<TH1F*>& enu_syst) {
  if (enu_syst.empty()) {
    std::cout << "EnuCollapsed: No multisims!?" << std::endl;
    return new TGraphErrors();
  }

  const size_t nbins = enu_syst[0]->GetNbinsX();

  // Compute the mean and standard deviation across universes using
  // Welford's method, cf. Art of Computer Programming (D. Knuth)
  double* x  = new double[nbins];
  double* y  = new double[nbins];
  double* y0 = new double[nbins];
  double* s  = new double[nbins];
  double* s0 = new double[nbins];
  double* xe = new double[nbins];

  for (size_t i=0; i<nbins; i++) {
    x[i] = enu_syst[0]->GetBinCenter(i);
    y[i] = y0[i] = enu_syst[0]->GetBinContent(i+1);
    s[i] = s0[i] = 0.0;

    for (size_t j=1; j<enu_syst.size(); j++) {
      double v = enu_syst[j]->GetBinContent(i);
      y[i] = y0[i] + (v - y0[i]) / j;
      s[i] = s0[i] + (v - y0[i]) * (v - y[i]);

      y0[i] = y[i];
      s0[i] = s[i];
    }

    xe[i] = 0.5 * enu_syst[0]->GetBinWidth(i+1);
    s[i] = sqrt(s[i] / enu_syst.size());
    std::cout << i << ": " << s[i] << std::endl;
  }

  TGraphErrors* g = new TGraphErrors(nbins, x, y, xe, s);

  delete[] x;
  delete[] y;
  delete[] y0;
  delete[] s;
  delete[] s0;
  delete[] xe;

  return g;
}

int WriteSystPlots(char* filename) {
  gStyle->SetOptStat(0);
  TFile* fin = NULL;

  //if (argc == 2) {
  //  fin = TFile::Open(argv[1]);
  //}
  //else {
  //  std::cerr << "Usage: " << argv[0] << " [input.root]" << std::endl;
  //  return 1;
  //}
  fin = TFile::Open(filename);
  assert(fin && fin->IsOpen());

  // Set up the input tree for reading
  TTree* _tree = (TTree*) fin->Get("data");
  assert(_tree && _tree->GetEntries() > 0);
  _tree->Print();

  galleryfmwk::TSSelection::OutputData _data;
  _tree->SetBranchAddress("enu", &_data.enu);
  _tree->SetBranchAddress("q2", &_data.q2);
  _tree->SetBranchAddress("w", &_data.w);
  _tree->SetBranchAddress("int", &_data.int_type);
  _tree->SetBranchAddress("mode", &_data.int_mode);
  _tree->SetBranchAddress("ccnc", &_data.ccnc);
  _tree->SetBranchAddress("eccqe", &_data.eccqe);
  _tree->SetBranchAddress("ep", &_data.ep);
  _tree->SetBranchAddress("ppdg", &_data.ppdg);
  _tree->SetBranchAddress("elep", &_data.elep);
  _tree->SetBranchAddress("lpdg", &_data.lpdg);
  _tree->SetBranchAddress("lpid", &_data.lpid);
  _tree->SetBranchAddress("bnbweight", &_data.bnbweight);
  _tree->SetBranchAddress("dataset", &_data.dataset);
  _tree->SetBranchAddress("weights", &_data.weights);

  // Histograms for universes
  std::map<std::string, std::vector<TH1F*> > hists;

  // Output file
  TFile f("fhist.root", "recreate");

  // Event loop
  long nentries = _tree->GetEntriesFast();
  size_t wnmin = 1e6;
  std::cout << "ENTRIES " << nentries << std::endl;

  nentries = 10000;
  bool init = false;
  for (long i=0; i<nentries; i++) {
    _tree->GetEntry(i);

    if (i % 1000 == 0) {
      std::cout << ".";
      std::cout.flush();
    }

    if (_data.lpid != 13 || _data.int_type != 1001) {
      continue;
    }
    
    // Create histograms for each weighting function and universe
    if (!init) {
      for (auto it=_data.weights->begin(); it!=_data.weights->end(); ++it) {
        std::string wname = it->first;
        //if (wname.find("genie") == std::string::npos && wname.find("xsr") == std::string::npos) {
        //  continue;
        //}
        if (wname.find("bnbcorrection") != std::string::npos) {
          continue;
        }
        size_t wn = it->second.size();
        hists[wname].resize(wn, NULL);
        std::cout << "resize " << wname << " -> " << wn << std::endl;
        for (size_t j=0; j<wn; j++) {
          char hname[50];
          snprintf(hname, 50, "h_%s_%03lu", wname.c_str(), j);
          std::cout << hname <<std::endl;
          //hists[wname][j] = new TH1F(hname, ";#nu Energy (GeV);Events", 15, 0, 3.0);
          hists[wname][j] = new TH1F(hname, ";CCQE Energy (MeV);Events", 25, 0, 2500);
        }

        if (wn < wnmin) {
            wnmin = wn;
        }
      }

      std::cout << "min " << wnmin <<std::endl;

      // Product
      std::string wname = "total";
      hists[wname].resize(wnmin, NULL);
      std::cout << "resize " << wname << " -> " << wnmin << std::endl;
      for (size_t j=0; j<wnmin; j++) {
        char hname[50];
        snprintf(hname, 50, "h_%s_%03lu", wname.c_str(), j);
        std::cout << hname <<std::endl;
        //hists[wname][j] = new TH1F(hname, ";#nu Energy (GeV);Events", 15, 0, 3.0);
        hists[wname][j] = new TH1F(hname, ";CCQE Energy (MeV);Events", 25, 0, 2500);
      }

      init = true;
    }

    std::vector<double> product(wnmin, 1.0);
    for (auto it=_data.weights->begin(); it!=_data.weights->end(); ++it) {
      std::string wname = it->first;
      size_t wn = it->second.size();
      if (wname.find("genie") == std::string::npos && wname.find("xsr") == std::string::npos) {
        continue;
      }
      //if (wname.find("qema") == std::string::npos) {
      //  continue;
      //}

      // Fill histograms
      //std::cout << hists[wname].size() << " " << it->second.size() << std::endl;
      for (size_t j=0; j<wn; j++) {
        //std::cout << wname << " " << j << " " << enu << " " << it->second[j] << std::endl;
        //hists[wname][j]->Fill(enu, it->second[j]);
        hists[wname][j]->Fill(_data.eccqe, it->second[j]);
        if (j < wnmin) {
          product[j] *= it->second[j];
        }
      }
    }

    for (size_t j=0; j<wnmin; j++) {
      //hists["total"][j]->Fill(enu, product[j]);
      hists["total"][j]->Fill(_data.eccqe, product[j]);
    }
  }

  // Write to file
  f.cd();
  for (auto it=hists.begin(); it!=hists.end(); ++it) {
    //if (it->first.find("qema") == std::string::npos) {
    //  continue;
    //}

    //for (size_t j=0; j<it->second.size(); j++) {
    //  hists[it->first][j]->Write();
    //}

    char hname[50];
    snprintf(hname, 50, "g_%s", it->first.c_str());
    TGraphErrors* g = EnuCollapsed2(it->second);
    g->SetName(hname);
    g->GetXaxis()->SetRangeUser(0, 3.0);
    g->SetTitle(";CCQE Energy (MeV);Events");
    g->SetFillColor(kRed);
    g->SetFillStyle(3001);
    g->Write();

    char cname[50];
    snprintf(cname, 50, "c_%s", it->first.c_str());
    TCanvas* c = new TCanvas(cname, "");
    c->cd();
    g->Draw("a2");
    //g->Draw("p");

    TPaveText* label = new TPaveText(1.5, 40, 3.0, 42.5);
    label->AddText(it->first.c_str());
    label->SetBorderSize(0);
    label->SetFillColor(kWhite);
    label->SetTextFont(132);
    label->Draw();

    c->Update();
    c->Write();
  }

  f.Close();
  fin->Close();

  return 0;
}

}  // namespace galleryfmwk

