/**
 * \file sel.h
 * \brief A truth-based event selection
 * \author A. Mastbaum <mastbaum@uchicago.edu>
 */

#ifndef GALLERY_FMWK_SEL_ANA_H
#define GALLERY_FMWK_SEL_ANA_H

#include <map>
#include <string>

#include "Analysis/ana_base.h"

class TFile;
class TH1F;
class TH2F;
class TNtuple;

namespace galleryfmwk {

/**
 * \class sel
 * \brief Truth-based selection approximating 1l1p
 */
class sel : galleryfmwk::ana_base {

public:
  sel() : _verbose(false) {}

  bool initialize();

  bool analyze(gallery::Event* ev);

  bool finalize();

  void setVerbose(bool b) { _verbose = b; }
  void setEventWeightProducer(std::string s) { _ew_producer = s; }
  void setMCTruthProducer(std::string s) { _mct_producer = s; }
  void setMCFluxProducer(std::string s) { _mcf_producer = s; }
  void setMCShowerProducer(std::string s) { _mcshw_producer = s; }
  void setMCTrackProducer(std::string s) { _mctrk_producer = s; }
  void setDatasetID(int id) { _dataset_id = id; }
    
protected:
  std::string _track_producer;
  std::string _ew_producer;
  std::string _mct_producer;
  std::string _mcf_producer;
  std::string _mctrk_producer;
  std::string _mcshw_producer;
  bool _verbose;
  int _dataset_id;

  TFile* _pdf_file;

  TNtuple* _data;

  size_t true_1e1p;
  size_t good_1e1p;
  size_t miss_1e1p;
  size_t true_1m1p;
  size_t good_1m1p;
  size_t miss_1m1p;

  std::map<int, TH2F*> _trackdedxs;
};

}  // namespace galleryfmwk

#endif  // GALLERY_FMWK_SEL_H

