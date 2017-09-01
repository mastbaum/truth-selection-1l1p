/**
 * \file TSPDFGen.h
 * \brief A truth-based event selection
 * \author A. Mastbaum <mastbaum@uchicago.edu>
 */

#ifndef GALLERY_FMWK_TSPDFGEN_H
#define GALLERY_FMWK_TSPDFGEN_H

#include <map>
#include <string>

#include "Analysis/ana_base.h"

class TH1F;
class TH2F;

namespace galleryfmwk {

/**
 * \class TSPDFGen
 * \brief Truth-based selection approximating 1l1p
 */
class TSPDFGen : galleryfmwk::ana_base {

public:
  TSPDFGen() : _verbose(false) {}

  bool initialize();

  bool analyze(gallery::Event* ev);

  bool finalize();

  void setVerbose(bool b) { _verbose = b; }
  void setEventWeightProducer(std::string s) { _ew_producer = s; }
  void setMCTruthProducer(std::string s) { _mct_producer = s; }
  void setMCFluxProducer(std::string s) { _mcf_producer = s; }
  void setMCShowerProducer(std::string s) { _mcshw_producer = s; }
  void setMCTrackProducer(std::string s) { _mctrk_producer = s; }
    
protected:
  std::string _track_producer;
  std::string _ew_producer;
  std::string _mct_producer;
  std::string _mcf_producer;
  std::string _mctrk_producer;
  std::string _mcshw_producer;
  bool _verbose;

  std::map<int, TH2F*> _trackdedxs;
  std::map<int, TH1F*> _showerdedxs;
};

}  // namespace galleryfmwk

#endif  // GALLERY_FMWK_TSPDFGEN_H

