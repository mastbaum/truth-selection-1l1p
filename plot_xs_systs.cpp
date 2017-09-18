/**
 * Plot distributions with systematics errors from weights.
 *
 * A. Mastbaum <mastbaum@uchicago.edu>, 2017/06/27
 */

#include <cassert>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>

struct OutputData {
  OutputData() {
    weights = NULL;
  }
  double enu;
  double q2;
  double w;
  double q0;
  double q3;
  int int_type;
  int int_mode;
  bool ccnc;
  double eccqe;
  double ep;
  int ppdg;
  double elep;
  int lpdg;
  int lpid;
  double bnbweight;
  int dataset;
  std::map<std::string, std::vector<double> >* weights;
};

// fixme
TGraphErrors* EnuCollapsed(std::vector<TH1F*>& enu_syst) {
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
    std::cout << i << ": " << y[i] << " " << s[i] << std::endl;
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

int main(int argc, char* argv[]) {
  gStyle->SetOptStat(0);
  TFile* fin = NULL;

  if (argc < 2) {
    fin = TFile::Open("/uboone/data/users/mastbaum/rw_test/fittree_bnb_rw2.root");
  }
  else if (argc == 2) {
    fin = TFile::Open(argv[1]);
  }
  else {
    std::cerr << "Usage: " << argv[0] << " [input.root]" << std::endl;
    return 1;
  }
  assert(fin && fin->IsOpen());

  // Set up the input tree for reading
  //TTree* ftree = (TTree*) fin->Get("fittree");
  //assert(ftree && ftree->GetEntries());

  //ftree->Print();

  //int nuPDG, intType;
  //double nuE, lepE;
  //std::map<std::string, std::vector<double> >* wgh = new std::map<std::string, std::vector<double> >;

  //ftree->SetBranchAddress("nuPDG", &nuPDG);
  //ftree->SetBranchAddress("type", &intType);
  //ftree->SetBranchAddress("nuEnergy", &nuE);
  //ftree->SetBranchAddress("leptonEnergy", &lepE);
  //ftree->SetBranchAddress("mcweight", &wgh);

  TTree* ftree = (TTree*) fin->Get("data");
  assert(ftree && ftree->GetEntries() > 0);
  ftree->Print();

  OutputData _data;
  ftree->SetBranchAddress("enu", &_data.enu);
  ftree->SetBranchAddress("q2", &_data.q2);
  ftree->SetBranchAddress("w", &_data.w);
  ftree->SetBranchAddress("int", &_data.int_type);
  ftree->SetBranchAddress("mode", &_data.int_mode);
  ftree->SetBranchAddress("ccnc", &_data.ccnc);
  ftree->SetBranchAddress("eccqe", &_data.eccqe);
  ftree->SetBranchAddress("ep", &_data.ep);
  ftree->SetBranchAddress("ppdg", &_data.ppdg);
  ftree->SetBranchAddress("elep", &_data.elep);
  ftree->SetBranchAddress("lpdg", &_data.lpdg);
  ftree->SetBranchAddress("lpid", &_data.lpid);
  ftree->SetBranchAddress("bnbweight", &_data.bnbweight);
  ftree->SetBranchAddress("dataset", &_data.dataset);
  ftree->SetBranchAddress("weights", &_data.weights);

  // Histograms for universes
  std::map<std::string, std::vector<TH1F*> > hists;

  // Output file
  TFile f("fhist.root", "recreate");

  // Event loop
  long nentries = ftree->GetEntriesFast();
  size_t wnmin = 1e6;
  std::cout << "ENTRIES " << nentries << std::endl;

  nentries = 10000;
  bool init = false;
  for (long i=0; i<nentries; i++) {
    ftree->GetEntry(i);

    if (i % 1000 == 0) {
      std::cout << ".";
      std::cout.flush();
    }

    if (_data.lpid != 13) {
      //std::cout << "got a pdg = " << nuPDG << ", type = " << intType << std::endl;
      continue;
    }
    
    // Create histograms for each weighting function and universe
    if (!init) {
      for (auto it=_data.weights->begin(); it!=_data.weights->end(); ++it) {
        std::string wname = it->first;
        if (wname.find("genie") == std::string::npos && wname.find("xsr") == std::string::npos) {
          continue;
        }
        if (wname.find("qema") == std::string::npos) {
          continue;
        }
        size_t wn = it->second.size();
        hists[wname].resize(wn, NULL);
        std::cout << "resize " << wname << " -> " << wn << std::endl;
        for (size_t j=0; j<wn; j++) {
          char hname[50];
          snprintf(hname, 50, "h_%s_%03i", wname.c_str(), j);
          std::cout << hname <<std::endl;
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
        snprintf(hname, 50, "h_%s_%03i", wname.c_str(), j);
        std::cout << hname <<std::endl;
        hists[wname][j] = new TH1F(hname, ";CCQE Energy (MeV);Events", 25, 0, 2500);
      }

      init = true;
    }

    std::vector<double> product(wnmin, 1.0);
    for (auto it=_data.weights->begin(); it!=_data.weights->end(); ++it) {
      std::string wname = it->first;
      size_t wn = it->second.size();
      //if (wname.find("genie") == std::string::npos && wname.find("xsr") == std::string::npos) {
      //  continue;
      //}
      if (wname.find("pionqex") == std::string::npos) {
        continue;
      }
      if (wname.find("bnb") != std::string::npos) {
        continue;
      }

      // Fill histograms
      //std::cout << hists[wname].size() << " " << it->second.size() << std::endl;
      for (size_t j=0; j<wn; j++) {
        if (hists.find(wname) == hists.end()) continue;
        if (it->second[j] > 3) {
          std::cout << wname << " " << j << " " << _data.enu << " " << it->second[j] << std::endl;
        }
        hists[wname][j]->Fill(_data.eccqe, it->second[j]);
        if (j < wnmin) {
          product[j] *= it->second[j];
        }
      }
    }

    for (size_t j=0; j<wnmin; j++) {
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
    TGraphErrors* g = EnuCollapsed(it->second);
    g->SetName(hname);
    g->GetXaxis()->SetRangeUser(0, 2500);
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

    TPaveText* label = new TPaveText(1500, 40, 2500, 42.5);
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

