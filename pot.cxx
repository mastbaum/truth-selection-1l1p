#include <iostream>

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "TFile.h"
#include "TTree.h"

#include "Analysis/ana_base.h"
#include "gallery/Event.h"
#include "pot.h"

namespace galleryfmwk {

bool pot::initialize() {
  totpot = 0;
  return true;
}

bool pot::analyze(gallery::Event* ev) {
  art::InputTag potsum_tag("generator");
  auto const& potsum = (*ev->getValidHandle<sumdata::POTSummary>(potsum_tag));
  totpot += potsum.totgoodpot;

  return true;
}


void pot::processFile(char* filename) {
  TFile f(filename);
  TTree* t = (TTree*) f.Get("SubRuns");
  Double_t pot;
  t->SetBranchAddress("sumdata::POTSummary_generator__GenieGen.obj.totgoodpot", &pot);
  t->GetEntry(0);
  totpot+=pot;
  std::cout << pot << " " << totpot << " " << filename << std::endl;
}


bool pot::finalize() {
  std::cout << "TOTAL POT: " << totpot << std::endl;

  return true;
}

}  // namespace gallery_fmwk

