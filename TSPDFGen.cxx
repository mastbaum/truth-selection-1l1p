#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "gallery/Event.h"

#include <TH2F.h>
#include <TH1F.h>

#include "TSUtil.h"
#include "TSPDFGen.h"
#include "TSSelection.h"

namespace galleryfmwk {

bool TSPDFGen::initialize() {
  return true;
}


bool TSPDFGen::analyze(gallery::Event* ev) {
  // Get handles for event data
  art::InputTag mctruth_tag(_mct_producer);
  auto const& mctruth_list = (*ev->getValidHandle<std::vector<simb::MCTruth> >(mctruth_tag));

  art::InputTag mcshower_tag(_mcshw_producer);
  auto const& mcshower_list = (*ev->getValidHandle<std::vector<sim::MCShower> >(mcshower_tag));

  art::InputTag mctrack_tag(_mctrk_producer);
  auto const& mctrack_list = (*ev->getValidHandle<std::vector<sim::MCTrack> >(mctrack_tag));

  _fout->cd();

  // Loop over MCTruth interactions
  for(size_t i=0; i<mctruth_list.size(); i++) {
    auto const& mctruth = mctruth_list.at(i);

    // Get vertex-associated contained tracks
    for (size_t j=0; j<mctrack_list.size(); j++) {
      const sim::MCTrack& mct = mctrack_list.at(j);

      if (!TSSelection::goodTrack(mct, mctruth)) {
        continue;
      }

      int pdg = mct.PdgCode();

      // Create dE/dx distribution for this PDG code
      if (_trackdedxs.find(pdg) == _trackdedxs.end()) {
        char hname[50];
        snprintf(hname, 50, "htrackdedx_%i", pdg);
        _trackdedxs[pdg] = new TH2F(hname, ";Residual range;dE/dx (MeV/cm?)", 100, 0, 200, 100, 0, 10);
      }

      assert(!mct.dEdx().empty());

      // Build up the track distribution
      TH2F* h = tsutil::HistDEdx(mct);
      _trackdedxs[pdg]->Add(h);
    }

    // Get vertex-associated contained showers
    for (size_t j=0; j<mcshower_list.size(); j++) {
      const sim::MCShower& mcs = mcshower_list.at(j);

      if (!TSSelection::goodShower(mcs, mctruth)) {
        continue;
      }

      int pdg = mcs.PdgCode();

      // Create dE/dx distribution
      if (_showerdedxs.find(pdg) == _showerdedxs.end()) {
        char hname[50];
        snprintf(hname, 50, "hshowerdedx_%i", pdg);
        _showerdedxs[pdg] = new TH1F(hname, ";dE/dx (MeV/cm?);Entries", 100, 0, 10);
      }

      // Build up the shower distribution
      double dedx = mcs.dEdx();
      _showerdedxs.at(mcs.PdgCode())->Fill(dedx);
    }
  }

  return true;
}


bool TSPDFGen::finalize() {
  // Write ROOT file output
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

