#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

struct OutputData {
  OutputData() {
    weights = NULL;
  }
  double enu;
  double q2;
  double w;
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

class CovarianceMatrix {
public:
  CovarianceMatrix();
  virtual ~CovarianceMatrix();

  void init();

  void analyze();

  /**
   * Container for an event sample, e.g. nue/numu/etc.
   */
  class EventSample {
    public:
      /**
       * Constructor.
       *
       * Note: You can optionally specify the number of systematics
       * universes up front with the nweights parameter. If this isn't
       * known until runtime, call Resize() later.
       *
       * \param _name String name for the event sample
       * \param nbins Number of energy bins
       * \param elo Lower limit of energy spectrum
       * \param ehi Upper limit of energy spectrum
       * \param nweights Number of systematics universes (see note)
       */
      EventSample(std::string _name="sample", size_t nbins=14,
                  double elo=0.2, double ehi=3.0, size_t nweights=0);

      /** Destructor. */
      ~EventSample();

      /**
       * Get the energy spectrum as a graph with error bars representing
       * the systematic uncertainty.
       */
      TGraphErrors* EnuCollapsed();

      /** Set the number of universes. */
      void Resize(size_t nweights);

      /**
       * Covariance Matrix
       *
       *   i && j     = energy bins
       *   n          = number of weights
       *   N^cv_i     = number of events in bin i for the central value
       *   N^syst_i,m = number of events in bin i for the systematic
       *                variation in universe "m"
       *   E_ij       = the covariance (square of the uncertainty) for
       *                bins i,j
       *   E_ij = (1/n) Sum((N^cv_i - N^syst_i,m)*(N^cv_j - N^syst_j,m), m)
       */
      static TH2D* CovarianceMatrix(TH1D* nom, std::vector<TH1D*> syst);

      /** Covariance matrix using internal histograms */
      TH2D* CovarianceMatrix();

      /** Correlation matrix: Corr[ij] = Cov[ij]/Sqrt(Cov[ii]*Cov[jj]) */
      static TH2D* CorrelationMatrix(TH2D* _cov);

      /** Correlation matrix using internal histograms */
      TH2D* CorrelationMatrix();

      std::string name;  //!< String name for this event sample
      TH1D* enu;  //!< "Nominal" energy spectrum
      std::vector<TH1D*> enu_syst;  //!< Spectra for each systematic universe

    protected:
      TH2D* cov;  //!< Cached covariance matrix
  };

private:
  std::set<std::string> use_weights;  //!< Weight functions to use
  std::vector<EventSample*> samples;  //!< Event samples
  TFile* fFile;  //!< File for output

  std::ofstream ofs;
};


/******************************************************************************
 ** CovarianceMatrix::EventSample implementation                             **
 *****************************************************************************/

CovarianceMatrix::EventSample::EventSample(std::string _name,
                                            size_t nbins,
                                            double elo, double ehi,
                                            size_t nweights)
    : name(_name), enu(nullptr), cov(nullptr) {
  enu = new TH1D(("enu_" + name).c_str(),
                 ";E_{#nu} [GeV];Entries per bin",
                 nbins, elo, ehi);

  Resize(nweights);
}


CovarianceMatrix::EventSample::~EventSample() {
  delete cov;
  delete enu;
}


TGraphErrors* CovarianceMatrix::EventSample::EnuCollapsed() {
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
    std::cout << "CovarianceMatrix::EventSample::EnuCollapsed: "
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


void CovarianceMatrix::EventSample::Resize(size_t nweights) {
  enu_syst.clear();
  for (size_t i=0; i<nweights; i++) {
    std::string hname = Form("enu_%s_%zu", name.c_str(), i);
    TH1D* h = (TH1D*) enu->Clone(hname.c_str());
    enu_syst.push_back(h);
  }
}


TH2D* CovarianceMatrix::EventSample::CovarianceMatrix(
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


TH2D* CovarianceMatrix::EventSample::CovarianceMatrix() {
  delete cov;
  cov = CovarianceMatrix(enu, enu_syst);
  cov->SetName(("cov_" + name).c_str());
  return cov;
}



TH2D* CovarianceMatrix::EventSample::CorrelationMatrix(TH2D* _cov) {
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


TH2D* CovarianceMatrix::EventSample::CorrelationMatrix() {
  // Compute the covariance matrix first, if we haven't already
  if (!cov) {
    CovarianceMatrix();
  }

  TH2D* cor = CorrelationMatrix(cov);
  cor->SetName(("cor_" + name).c_str());
  return cor;
}


/******************************************************************************
 ** CovarianceMatrix implementation                                          **
 *****************************************************************************/

CovarianceMatrix::CovarianceMatrix() {
  fFile = TFile::Open("cov.root", "recreate");
  assert(fFile);

  init();
}


CovarianceMatrix::~CovarianceMatrix() {
  ofs.close();
}


void CovarianceMatrix::init() {
  samples.push_back(new EventSample("numu"));
  samples.push_back(new EventSample("nue"));

  std::vector<std::string> wv = {
    // "*"
    //,"bnbcorrection_FluxHist"
     "expskin_FluxUnisim"
    ,"genie_AGKYpT_Genie"
    ,"genie_AGKYxF_Genie"
    ,"genie_DISAth_Genie"
    ,"genie_DISBth_Genie"
    ,"genie_DISCv1u_Genie"
    ,"genie_DISCv2u_Genie"
    ,"genie_FermiGasModelKf_Genie"
    ,"genie_FermiGasModelSf_Genie"
    ,"genie_FormZone_Genie"
    ,"genie_IntraNukeNabs_Genie"
    ,"genie_IntraNukeNcex_Genie"
    ,"genie_IntraNukeNel_Genie"
    ,"genie_IntraNukeNinel_Genie"
    ,"genie_IntraNukeNmfp_Genie"
    ,"genie_IntraNukeNpi_Genie"
    ,"genie_IntraNukePIabs_Genie"
    ,"genie_IntraNukePIcex_Genie"
    ,"genie_IntraNukePIel_Genie"
    ,"genie_IntraNukePIinel_Genie"
    ,"genie_IntraNukePImfp_Genie"
    ,"genie_IntraNukePIpi_Genie"
    ,"genie_NC_Genie"
    ,"genie_NonResRvbarp1pi_Genie"
    ,"genie_NonResRvbarp2pi_Genie"
    ,"genie_NonResRvp1pi_Genie"
    ,"genie_NonResRvp2pi_Genie"
    ,"genie_ResDecayEta_Genie"
    ,"genie_ResDecayGamma_Genie"
    ,"genie_ResDecayTheta_Genie"
    ,"genie_ccresAxial_Genie"
    ,"genie_ccresVector_Genie"
    ,"genie_cohMA_Genie"
    ,"genie_cohR0_Genie"
    ,"genie_ncelAxial_Genie"
    ,"genie_ncelEta_Genie"
    ,"genie_ncresAxial_Genie"
    ,"genie_ncresVector_Genie"
    ,"genie_qema_Genie"
    ,"genie_qevec_Genie"
    ,"horncurrent_FluxUnisim"
    ,"kminus_PrimaryHadronNormalization"
    ,"kplus_PrimaryHadronFeynmanScaling"
    ,"kzero_PrimaryHadronSanfordWang"
    ,"nucleoninexsec_FluxUnisim"
    ,"nucleonqexsec_FluxUnisim"
    ,"nucleontotxsec_FluxUnisim"
    ,"piminus_PrimaryHadronSWCentralSplineVariation"
    ,"pioninexsec_FluxUnisim"
    ,"pionqexsec_FluxUnisim"
    ,"piontotxsec_FluxUnisim"
    ,"piplus_PrimaryHadronSWCentralSplineVariation"
  };

  // "Fv3technote_XSecRatio", "*"

  use_weights = std::set<std::string>(wv.begin(), wv.end());

  std::string out_file = "cov.txt";
  ofs.open(out_file, std::ofstream::out);

  std::cout << "CovarianceMatrix: Initialized. Weights: ";
  for (auto it : use_weights) {
    std::cout << it << " ";
  }
  std::cout << std::endl;

  std::cout << "CovarianceMatrix: Writing output to " << out_file << std::endl;
}
  

void CovarianceMatrix::analyze() {
  // Grab the MCTruth information
  //TFile f("sel/t_e.root");
  TFile f("./sel_num_nue.root");
  TTree* _tree = (TTree*) f.Get("data");
  assert(_tree && _tree->GetEntries() > 0);

  OutputData _data;
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

      std::cout << "CovarianceMatrix: Found weight: " << it.first << std::endl;

      // Compute universe-wise product of all requsted weights
      if (use_weights.find("*") != use_weights.end() ||
          use_weights.find(it.first) != use_weights.end()) {
        for (size_t i=0; i<weights.size(); i++) {
          weights[i] *= it.second[i];
        }
      }
    }

    // Neutrino interaction truth
    double nuEnergy = _data.enu;

    // Determine which event sample this event corresponds to
    EventSample* sample = nullptr;

    if (_data.lpdg == 11) {
      sample = samples[1];
    }
    else if (_data.lpdg == 13) {
      sample = samples[0];
    }
    else {
      std::cout << "Unknown lepton PDG " << _data.lpdg << std::endl;
      continue;
    }

    // Fill histograms for this event sample
    if (sample->enu_syst.empty()) {
      sample->Resize(weights.size());
    }
    else {
      assert(sample->enu_syst.size() == weights.size());
    }

    // Neutrino Energy
    sample->enu->Fill(nuEnergy);

    // Fill weights
    std::cout << "CovarianceMatrix: Weights: ";
    ofs << _data.lpdg << "\t" << nuEnergy << "\t";
    for (size_t i=0; i<weights.size(); i++) {
      sample->enu_syst[i]->Fill(nuEnergy, weights[i]);
      std::cout << weights[i] << " ";
      ofs << weights[i] << "\t";
    }
    std::cout << std::endl;
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


int main(int argc, char* argv[]) {
  CovarianceMatrix cov;
  cov.analyze();
  return 0;
}

