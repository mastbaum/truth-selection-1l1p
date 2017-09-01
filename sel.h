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

class TDatabasePDG;
class TFile;
class TH2F;
class TNtuple;
class TTree;

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

  // Set the producers for data products
  void setFluxWeightProducer(std::string s) { _fw_producer = s; }
  void setEventWeightProducer(std::string s) { _ew_producer = s; }
  void setMCTruthProducer(std::string s) { _mct_producer = s; }
  void setMCFluxProducer(std::string s) { _mcf_producer = s; }
  void setMCShowerProducer(std::string s) { _mcshw_producer = s; }
  void setMCTrackProducer(std::string s) { _mctrk_producer = s; }

  // Set a numeric dataset ID, which is written into the tree
  void setDatasetID(int id) { _dataset_id = id; }

  double get_mass(int pdg) const;

  struct OutputData {
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
    std::map<std::string, std::vector<double> > weights;
  };
    
protected:
  // Data product producers
  std::string _track_producer;
  std::string _fw_producer;
  std::string _ew_producer;
  std::string _mct_producer;
  std::string _mcf_producer;
  std::string _mctrk_producer;
  std::string _mcshw_producer;

  // Counters for efficiency and purity calculations
  size_t true_1e1p;
  size_t good_1e1p;
  size_t miss_1e1p;
  size_t true_1m1p;
  size_t good_1m1p;
  size_t miss_1m1p;

  bool _verbose;  //!< Print verbose output
  int _dataset_id;  //!< An arbitrary numeric ID
  OutputData* _data;  //!< Output data
  TTree* _tree;  //!< Output tree

  TNtuple* _truthtree;
  TNtuple* _mectree;
  TDatabasePDG* _pdgtable;

  TFile* _pdf_file;  //!< File containing dE/dx PDFs
  std::map<int, TH2F*> _trackdedxs;  // Track dE/dx distributions
};

}  // namespace galleryfmwk

#endif  // GALLERY_FMWK_SEL_H

