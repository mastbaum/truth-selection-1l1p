#include <cassert>
#include <iostream>
#include <sstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "gallery/Event.h"

#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TNtuple.h>

#include "tsutil.h"
#include "sel.h"

namespace galleryfmwk {

bool sel::initialize() {
  // Load track dE/dx distributions from file
  _pdf_file = TFile::Open("./pdfs.root");
  assert(_pdf_file->IsOpen());

  TIter next(_pdf_file->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    const char* name = key->GetName();
    TObjArray* tokens = TString(name).Tokenize("_");
    TString htype = ((TObjString*)tokens->At(0))->GetString();
    int pdg = atoi(((TObjString*)tokens->At(1))->GetString());

    if (htype.Contains("htrackdedx")) {
      TH2F* h = (TH2F*) _pdf_file->Get(name);
      if (h->Integral() == 0 || pdg < 0 || pdg > 10000) { continue; }
      std::cout << name << ": " << pdg << " " << key->GetClassName() << std::endl;

      // Ignore low bins
      for (int i=0; i<h->GetNbinsX(); i++) {
        for (int j=0; j<h->GetNbinsY(); j++) {
          if (i < 2 || j < 2) {
            h->SetBinContent(i, j, 0);
          }
        }
      }
      _trackdedxs[pdg] = h;
    }
  }

  _dataset_id = -1;
  true_1e1p = 0;
  good_1e1p = 0;
  miss_1e1p = 0;
  true_1m1p = 0;
  good_1m1p = 0;
  miss_1m1p = 0;

  if (_fout) {
    _fout->cd();
    _data = new TNtuple("data", "", "enu:int:mode:ccnc:q2:w:eccqe:ep:ppdg:elep:lpdg:lpid:dataset");
  }

  return true;
}

struct PIDParticle {
  int pdg;
  int pdgtrue;
  TLorentzVector p;
  double evis;
  double eccqe;

  friend std::ostream& operator<<(std::ostream& os, const PIDParticle& dt);
};

std::ostream& operator<<(std::ostream& os, const PIDParticle& dt) {
  os << dt.pdg << "(" << dt.evis << ")";
  return os;
}

bool is1l1p(std::vector<PIDParticle>& p, int lpdg) {
  if (p.size() > 2) {
    return false;
  }
  size_t np = 0;
  size_t nl = 0;
  for (size_t i=0; i<p.size(); i++) {
    if (p[i].pdg == 2212) {
      np++;
    }
    if (p[i].pdg == lpdg) {
      nl++;
    }
  }
  return np == 1 && nl == 1;
}

bool sel::analyze(gallery::Event* ev) {
  // Get event data
  // art::InputTag eventweight_tag(_ew_producer);
  // auto const& eventweights_list = (*ev->getValidHandle<std::vector<evwgh::MCEventWeight> >(eventweight_tag));

  art::InputTag mctruth_tag(_mct_producer);
  auto const& mctruth_list = (*ev->getValidHandle<std::vector<simb::MCTruth> >(mctruth_tag));

  art::InputTag mcshower_tag(_mcshw_producer);
  auto const& mcshower_list = (*ev->getValidHandle<std::vector<sim::MCShower> >(mcshower_tag));

  art::InputTag mctrack_tag(_mctrk_producer);
  auto const& mctrack_list = (*ev->getValidHandle<std::vector<sim::MCTrack> >(mctrack_tag));

  _fout->cd();

  for(size_t i=0; i<mctruth_list.size(); i++) {
    auto const& mctruth = mctruth_list.at(i);

    size_t ntracks = 0, nshowers = 0;

    // Keep track of particle content
    std::vector<PIDParticle> particles_found;
    std::vector<PIDParticle> particles_true;

    // Get vertex-associated contained tracks
    for (size_t j=0; j<mctrack_list.size(); j++) {
      const sim::MCTrack& mct = mctrack_list.at(j);
      if (mct.empty() || !isFromNuVertex(mctruth, mct) ||
          mct.Start().E() < 60 || !inFV(mct) ) {
        continue;
      }

      ntracks++;

      particles_true.push_back({
        mct.PdgCode(),
        mct.PdgCode(),
        mct.Start().Momentum(),
        mct.Start().E(),
        eccqe(mct.Start().Momentum())
      });

      // Build a distribution of dE/dx vs. residual range
      TH2F* htemp = new TH2F("htemp", ";Residual range (cm?);dE/dx (MeV/cm)",
                             100, 0, 200, 100, 0, 10);
      double s = 0;  // Track length
      TLorentzVector pos = mct.End().Position();  // Start from end, work back
      for (long k=mct.size()-1; k>=0; k--) {
        double dedx = mct.dEdx()[k];
        if (htemp->GetXaxis()->FindBin(s) >= 2 &&
            htemp->GetYaxis()->FindBin(dedx) >= 2) {
          htemp->Fill(s, dedx);
        }
        s += (pos.Vect() - mct[k].Position().Vect()).Mag();
        pos = mct[k].Position();
      }

      // KS test to pick the best match for PID
      double ks_best = 0;
      int pdg_best = -999;
      for (auto const& it : _trackdedxs) {
        if (htemp->Integral() == 0 || it.second->Integral() == 0) { continue; }
        double ks = it.second->KolmogorovTest(htemp);
        if (ks > ks_best) {
          ks_best = ks;
          pdg_best = it.first;
        }
      }

      delete htemp;

      // Un-PID "protons" that are too long or short
      if (pdg_best == 2212 && (s > 80 || s < 12)) {
        pdg_best = -888;
      } 

      // Call all remaining unmatched tracks protons
      if (pdg_best == -999) {
        pdg_best = 2212;
      }

      particles_found.push_back({
        pdg_best,
        mct.PdgCode(),
        mct.Start().Momentum(),
        mct.Start().E(),
        eccqe(mct.Start().Momentum())
      });
    }

    // Get vertex-associated contained showers
    for (size_t j=0; j<mcshower_list.size(); j++) {
      const sim::MCShower& mcs = mcshower_list.at(j);
      if (!isFromNuVertex(mctruth, mcs) ||
          mcs.Start().E() < 30 || !inFV(mcs)) {
        continue;
      }

      nshowers++;

      particles_true.push_back({
        mcs.PdgCode(),
        mcs.PdgCode(),
        mcs.Start().Momentum(),
        mcs.Start().E(),
        eccqe(mcs.Start().Momentum())
      });

      // Guess the PDG based on shower dE/dx (cut at 3.5)
      int pdg_best = (mcs.dEdx() < 3.5 ? 11 : 22);

      particles_found.push_back({
        pdg_best,
        mcs.PdgCode(),
        mcs.Start().Momentum(),
        mcs.Start().E(),
        eccqe(mcs.Start().Momentum())
      });
    }

    // Classify the event (found/true 1l1p/1m1p)
    bool f_1e1p = is1l1p(particles_found, 11);
    bool t_1e1p = is1l1p(particles_true, 11);
    bool f_1m1p = is1l1p(particles_found, 13);
    bool t_1m1p = is1l1p(particles_true, 13);

    if (t_1e1p) true_1e1p++;
    if (f_1e1p && t_1e1p) good_1e1p++;
    if (f_1e1p && !t_1e1p) miss_1e1p++;

    if (t_1m1p) true_1m1p++;
    if (f_1m1p && t_1m1p) good_1m1p++;
    if (f_1m1p && !t_1m1p) miss_1m1p++;

    // Print mis-IDs
    if ((f_1e1p && !t_1e1p) || (f_1m1p && !t_1m1p)) {
      std::cout << "true: [" << mctruth.GetNeutrino().InteractionType() << "] ";
      for (size_t k=0; k<particles_true.size(); k++) {
        std::cout << particles_true[k] << " ";
      }
      std::cout << std::endl;
      std::cout << "est: " << ntracks << " tracks, " << nshowers << " showers; ";
      for (size_t k=0; k<particles_found.size(); k++) {
        std::cout << particles_found[k] << " ";
      }
      std::cout << std::endl;
    }

    // Write event to output tree
    if (f_1e1p || f_1m1p) {
      double eccqe, ep, ppdg, elep, lpdg, lpid;

      for (size_t k=0; k<particles_found.size(); k++) {
        if (particles_found[k].pdg == 2212) {
          ep = particles_found[k].evis;
          ppdg = particles_found[k].pdgtrue;
        }
        else {
          eccqe = particles_found[k].eccqe;
          elep = particles_found[k].evis;
          lpid = particles_found[k].pdg;
          lpdg = particles_found[k].pdgtrue;
        }
      }

      const simb::MCNeutrino& nu = mctruth.GetNeutrino();
      _data->Fill(nu.Nu().E(), nu.InteractionType(), nu.Mode(), nu.CCNC(), nu.QSqr(), nu.W(),
                  eccqe, ep, ppdg, elep, lpdg, lpid, _dataset_id);
    }
  }

  return true;
}


bool sel::finalize() {
  std::cout << "1e1p true: " << true_1e1p << ", good: " << good_1e1p << ", miss: " << miss_1e1p << std::endl;
  std::cout << "1e1p eff: " << 1.0 * (good_1e1p + miss_1e1p) / true_1e1p << std::endl;
  std::cout << "1e1p pur: " << 1.0 * good_1e1p / (good_1e1p + miss_1e1p) << std::endl;

  std::cout << "1m1p true: " << true_1m1p << ", good: " << good_1m1p << ", miss: " << miss_1m1p << std::endl;
  std::cout << "1m1p eff: " << 1.0 * (good_1m1p + miss_1m1p) / true_1m1p << std::endl;
  std::cout << "1m1p pur: " << 1.0 * good_1m1p / (good_1m1p + miss_1m1p) << std::endl;

  if (_fout) {
    _fout->cd();
    _data->Write();
  }

  _pdf_file->Close();

  return true;
}

}  // namespace gallery_fmwk

