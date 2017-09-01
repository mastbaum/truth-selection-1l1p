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
#include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "uboone/EventWeight/MCEventWeight.h"
#include "gallery/Event.h"

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TTree.h>

#include "tsutil.h"
#include "sel.h"

namespace galleryfmwk {

// A struct to hold temporary track/shower data
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


// Utility function to test if a list of particles is 1l1p
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


// Based on https://cdcvs.fnal.gov/redmine/issues/15071
double sel::get_mass(int pdg) const {
  if (pdg < 1000000000) {
    TParticlePDG* ple = _pdgtable->GetParticle(pdg);
    return ple->Mass() * 1000.0;
  }
  else {
    int p = (pdg % 10000000) / 10000;
    int n = (pdg % 10000) / 10 - p;
    return (_pdgtable->GetParticle(2212)->Mass() * p + _pdgtable->GetParticle(2112)->Mass() * n) * 1000.0;
  }
}


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
      //std::cout << name << ": " << pdg << " " << key->GetClassName() << std::endl;

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

  // PDG table
  _pdgtable = new TDatabasePDG;

  // Initialize dataset identifier
  _dataset_id = -1;

  // Initialize event counters
  true_1e1p = 0;
  good_1e1p = 0;
  miss_1e1p = 0;
  true_1m1p = 0;
  good_1m1p = 0;
  miss_1m1p = 0;

  // Set up the output trees
  if (_fout) {
    _fout->cd();
    _data = new OutputData;
    _tree = new TTree("data", "");
    _tree->Branch("enu", &_data->enu);
    _tree->Branch("q2", &_data->q2);
    _tree->Branch("w", &_data->w);
    _tree->Branch("q0", &_data->q0);
    _tree->Branch("q3", &_data->q3);
    _tree->Branch("int", &_data->int_type);
    _tree->Branch("mode", &_data->int_mode);
    _tree->Branch("ccnc", &_data->ccnc);
    _tree->Branch("eccqe", &_data->eccqe);
    _tree->Branch("ep", &_data->ep);
    _tree->Branch("ppdg", &_data->ppdg);
    _tree->Branch("elep", &_data->elep);
    _tree->Branch("lpdg", &_data->lpdg);
    _tree->Branch("lpid", &_data->lpid);
    _tree->Branch("bnbweight", &_data->bnbweight);
    _tree->Branch("dataset", &_data->dataset);
    _tree->Branch("weights", &_data->weights);

    _truthtree = new TNtuple("truth", "", "nupdg:enu:ccnc:int:mode:w:q2:lpdg:elep:tlep:npip:npim:npi0:np:nn:fw:ttrk:rtrk:texit:tshr:rshr:sexit");
    _mectree = new TNtuple("mec", "", "nupdg:enu:ccnc:mode:w:q2:lpdg:tlep:ep0:ep1:ep2:ep3:ep4");
  }

  return true;
}


bool sel::analyze(gallery::Event* ev) {
  // Get event data
  art::InputTag gtruth_tag(_mct_producer);
  auto const& gtruth_list = (*ev->getValidHandle<std::vector<simb::GTruth> >(gtruth_tag));

  //art::InputTag fluxweight_tag(_fw_producer);
  //auto const& fluxweights_list = (*ev->getValidHandle<std::vector<evwgh::MCEventWeight> >(fluxweight_tag));

  art::InputTag eventweight_tag(_ew_producer);
  auto const& eventweights_list = (*ev->getValidHandle<std::vector<evwgh::MCEventWeight> >(eventweight_tag));

  art::InputTag mctruth_tag(_mct_producer);
  auto const& mctruth_list = (*ev->getValidHandle<std::vector<simb::MCTruth> >(mctruth_tag));

  art::InputTag mcshower_tag(_mcshw_producer);
  auto const& mcshower_list = (*ev->getValidHandle<std::vector<sim::MCShower> >(mcshower_tag));

  art::InputTag mctrack_tag(_mctrk_producer);
  auto const& mctrack_list = (*ev->getValidHandle<std::vector<sim::MCTrack> >(mctrack_tag));

  _fout->cd();

  assert(mctruth_list.size() == gtruth_list.size());
  assert(!eventweights_list.empty()); // && !fluxweights_list.empty());

  //for (size_t i=0; i<fluxweights_list.size(); i++) {
  //  std::cout << "Flux MCEventWeights: ";
  //  for (auto const& it : fluxweights_list[i].fWeight) {
  //    std::cout << it.first << "[" << it.second.size() << "] ";
  //  }
  //  std::cout << std::endl;
  //}

  //for (size_t i=0; i<eventweights_list.size(); i++) {
  //  std::cout << "Other MCEventWeights: ";
  //  for (auto const& it : eventweights_list[i].fWeight) {
  //    std::cout << it.first << "[" << it.second.size() << "] ";
  //  }
  //  std::cout << std::endl;
  //}

  for(size_t i=0; i<mctruth_list.size(); i++) {
    const simb::MCTruth& mctruth = mctruth_list.at(i);
    const simb::GTruth& gtruth = gtruth_list.at(i);  // Hope so...

    size_t ntracks = 0, nshowers = 0, ntexiting = 0, nsexiting = 0;

    // Keep track of particle content
    std::vector<PIDParticle> particles_found;
    std::vector<PIDParticle> particles_true;

    // Get vertex-associated contained tracks
    for (size_t j=0; j<mctrack_list.size(); j++) {
      const sim::MCTrack& mct = mctrack_list.at(j);
      if (!inFV(mct)) { ntexiting++; }

      //std::cout << mct.PdgCode() << " " << get_mass(mct.PdgCode()) << std::endl;
      if (mct.empty() || !isFromNuVertex(mctruth, mct) || !(mct.Process() == "primary")) {
        continue;
      }

      if (mct.Start().E() - get_mass(mct.PdgCode()) < 60) { // || !inFV(mct) ) {
        continue;
      }

      particles_true.push_back({
        mct.PdgCode(),
        mct.PdgCode(),
        mct.Start().Momentum(),
        mct.Start().E() - get_mass(mct.PdgCode()),
        eccqe(mct.Start().Momentum())
      });


      ntracks++;

      // Build a distribution of dE/dx vs. residual range
      TH2F* htemp = new TH2F("htemp", ";Residual range (cm?);dE/dx (MeV/cm)",
                             100, 0, 200, 100, 0, 10);
      double s = 0;  // Track length
      TLorentzVector pos = mct.End().Position();  // Start from end, work back
      for (long k=mct.size()-2; k>=0; k--) {
        //std::cout << mct.size() << " " << mct.dEdx().size() << " " << k << std::endl;
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
        mct.Start().E() - get_mass(mct.PdgCode()),
        eccqe(mct.Start().Momentum())
      });
    }

    // Get vertex-associated contained showers
    for (size_t j=0; j<mcshower_list.size(); j++) {
      const sim::MCShower& mcs = mcshower_list.at(j);
      if (!inFV(mcs)) { nsexiting++; }

      if (!isFromNuVertex(mctruth, mcs) || !(mcs.Process() == "primary")) {
        continue;
      }

      if (mcs.Start().E() - get_mass(mcs.PdgCode()) < 30) { // || !inFV(mcs)) {
        continue;
      }

      particles_true.push_back({
        mcs.PdgCode(),
        mcs.PdgCode(),
        mcs.Start().Momentum(),
        mcs.Start().E() - get_mass(mcs.PdgCode()),
        eccqe(mcs.Start().Momentum())
      });

      nshowers++;

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
    if ((f_1e1p && !t_1e1p) || (f_1m1p && !t_1m1p) /*|| mctruth.GetNeutrino().Mode() == 10*/) {
      std::cout << "true: " << mctruth.GetNeutrino().Nu().E() * 1000 << "[" << mctruth.GetNeutrino().InteractionType() << "] ";
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
      double eccqe=-1, ep=-1, ppdg=-1, elep=-1, lpdg=-1, lpid=-1;

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
      const simb::MCParticle& pnu = nu.Nu();
      const simb::MCParticle& plep = nu.Lepton();
      TLorentzVector xp = (pnu.Momentum() - plep.Momentum());
      
      _data->enu = nu.Nu().E();
      _data->q2 = nu.QSqr();
      _data->w = nu.W();
      _data->q0 = xp.E();
      _data->q3 = xp.Vect().Mag();
      _data->int_type = nu.InteractionType();
      _data->int_mode = nu.Mode();
      _data->ccnc = nu.CCNC();
      _data->eccqe = eccqe;
      _data->ep = ep;
      _data->ppdg = ppdg;
      _data->elep = elep;
      _data->lpdg = lpdg;
      _data->lpid = lpid;
      _data->bnbweight = 1.0; //fluxweights_list[0].fWeight.at("bnbcorrection_FluxHist").at(0);
      _data->dataset = _dataset_id;
      _data->weights = eventweights_list[0].fWeight;
      _tree->Fill();
    }

    // Fill the event truth tree
    float vtt[22] = {
      (float) mctruth.GetNeutrino().Nu().PdgCode(),
      (float) mctruth.GetNeutrino().Nu().E(),
      (float) mctruth.GetNeutrino().CCNC(),
      (float) mctruth.GetNeutrino().InteractionType(),
      (float) mctruth.GetNeutrino().Mode(),
      (float) mctruth.GetNeutrino().W(),
      (float) mctruth.GetNeutrino().QSqr(),
      (float) mctruth.GetNeutrino().Lepton().PdgCode(),
      (float) mctruth.GetNeutrino().Lepton().E(),
      (float) gtruth.fGint,
      (float) gtruth.fNumPiPlus,
      (float) gtruth.fNumPiMinus,
      (float) gtruth.fNumPi0,
      (float) gtruth.fNumProton,
      (float) gtruth.fNumNeutron,
      (float) 1.0, //fluxweights_list[0].fWeight.at("bnbcorrection_FluxHist").at(0),
      (float) mctrack_list.size(),
      (float) ntracks,
      (float) ntexiting,
      (float) mcshower_list.size(),
      (float) nshowers,
      (float) nsexiting
    };
    _truthtree->Fill(vtt);

    // MEC proton energies
    if (mctruth.GetNeutrino().Mode() == simb::kMEC) {
      //double ketot = 0;
      //std::cout << "MEC" << std::endl;
      std::vector<float> epmec(5, -1);
      for (size_t k=0; k<mctrack_list.size(); k++) {
        const sim::MCTrack& t = mctrack_list[k];
        if (t.PdgCode() == 2212 && t.Process() == "primary") {
          double ke = t.Start().E() - get_mass(t.PdgCode());
          epmec.push_back(ke);
          //ketot += ke / 1000;
          //std::cout << "TRK " << k << " KE " << ke
          //          << " ID " << t.TrackID() << " PROC " << t.Process()
          //          << " MID " << t.MotherTrackID() << " MPROC " << t.MotherProcess()
          //          << " AID " << t.AncestorTrackID() << " APROC " << t.AncestorProcess()
          //          << std::endl;
        }
      }

      std::sort(epmec.begin(), epmec.end(), std::greater<>());

      float v2[13] = {
        (float) mctruth.GetNeutrino().Nu().PdgCode(),
        (float) mctruth.GetNeutrino().Nu().E(),
        (float) mctruth.GetNeutrino().CCNC(),
        (float) mctruth.GetNeutrino().Mode(),
        (float) mctruth.GetNeutrino().W(),
        (float) mctruth.GetNeutrino().QSqr(),
        (float) mctruth.GetNeutrino().Lepton().PdgCode(),
        (float) mctruth.GetNeutrino().Lepton().E(),
        (float) epmec[0],
        (float) epmec[1],
        (float) epmec[2],
        (float) epmec[3],
        (float) epmec[4]
      };
      //std::cout << "ENU " << mctruth.GetNeutrino().Nu().E() << " ETOT " << ketot
      //          << (ketot > mctruth.GetNeutrino().Nu().E() ? " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" : "")
      //          << std::endl;
      _mectree->Fill(v2);
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
    _tree->Write();
    _truthtree->Write();
    _mectree->Write();
  }

  return true;
}

}  // namespace gallery_fmwk

