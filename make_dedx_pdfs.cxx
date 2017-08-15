#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "gallery/Event.h"

#include <TLorentzVector.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TNtuple.h>

#include "tsutil.h"
#include "make_dedx_pdfs.h"

namespace galleryfmwk {

bool make_dedx_pdfs::initialize() {
  return true;
}


bool make_dedx_pdfs::analyze(gallery::Event* ev) {
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

    // Get vertex-associated contained tracks
    for (size_t j=0; j<mctrack_list.size(); j++) {
      const sim::MCTrack& mct = mctrack_list.at(j);
      if (mct.empty() || mct.Start().E() < 60 ||
          !isFromNuVertex(mctruth, mct) || !inFV(mct)) {
        continue;
      }

      int pdg = mct.PdgCode();

      if (_trackdedxs.find(pdg) == _trackdedxs.end()) {
        char hname[50];
        snprintf(hname, 50, "htrackdedx_%i", pdg);
        _trackdedxs[pdg] = new TH2F(hname, ";Residual range;dE/dx (MeV/cm?)", 100, 0, 200, 100, 0, 10);
      }

      assert(!mct.dEdx().empty());

      double s = 0;
      TLorentzVector pos = mct.End().Position();
      for (long k=mct.size()-1; k>=0; k--) {
        double dedx = mct.dEdx()[k];
        _trackdedxs.at(pdg)->Fill(s, dedx);
        s += (pos.Vect() - mct[k].Position().Vect()).Mag();
        pos = mct[k].Position();
      }
    }

    // Get vertex-associated contained showers
    for (size_t j=0; j<mcshower_list.size(); j++) {
      const sim::MCShower& mcs = mcshower_list.at(j);
      if (!isFromNuVertex(mctruth, mcs) || !inFV(mcs)) {
        continue;
      }

      int pdg = mcs.PdgCode();

      if (_showerdedxs.find(pdg) == _showerdedxs.end()) {
        char hname[50];
        snprintf(hname, 50, "hshowerdedx_%i", pdg);
        _showerdedxs[pdg] = new TH1F(hname, ";dE/dx (MeV/cm?);Entries", 100, 0, 10);
      }

      double dedx = mcs.dEdx();
      _showerdedxs.at(mcs.PdgCode())->Fill(dedx);
    }
  }

  return true;
}


bool make_dedx_pdfs::finalize() {
  if (_fout) {
    _fout->cd();

    for (auto& it : _trackdedxs) {
      it.second->Write();
    }

    for (auto& it : _showerdedxs) {
      it.second->Write();
    }
  }

  return true;
}

}  // namespace gallery_fmwk

