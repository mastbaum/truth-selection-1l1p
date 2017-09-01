#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "Analysis/ana_base.h"
#include "gallery/Event.h"

namespace galleryfmwk {

class pot : galleryfmwk::ana_base {
public:
  pot() {}
  bool initialize();
  bool analyze(gallery::Event* ev);
  bool finalize();

  void processFile(char* filename);

protected:
  double totpot;
};

}  // namespace gallery_fmwk

