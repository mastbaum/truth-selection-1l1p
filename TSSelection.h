/**
 * \file TSSelection.h
 * \brief A truth-based 1l1p event selection
 * \author A. Mastbaum <mastbaum@uchicago.edu>
 */

#ifndef GALLERY_FMWK_TSSELECTION_H
#define GALLERY_FMWK_TSSELECTION_H

#include <map>
#include <string>
#include <random>

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "gallery/Event.h"
#include "TSUtil.h"

class TDatabasePDG;
class TFile;
class TH2F;
class TNtuple;
class TTree;

namespace galleryfmwk {

/**
 * \class TSSelection
 * \brief Truth-based selection approximating 1l1p
 */
class TSSelection {
  enum EventType {
    P0, P1, PN, TRKN, ANY
  };

  struct EventCounts {
    size_t true_1e;
    size_t good_1e;
    size_t miss_1e;
    size_t true_1m;
    size_t good_1m;
    size_t miss_1m;


    inline void fill(bool f_1e, bool t_1e, bool f_1m, bool t_1m) {
      if (t_1e) true_1e ++;
      if (f_1e && t_1e) good_1e++;
      if (f_1e && !t_1e) miss_1e++;
      if (t_1m) true_1m++;
      if (f_1m && t_1m) good_1m++;
      if (f_1m && !t_1m) miss_1m++;
    }
  };

public:
  // A structure to hold temporary track/shower data during processing
  struct PIDParticle {
    int pdg;
    int pdgtrue;
    TLorentzVector p;
    double evis;
    double eccqe;
    double len;
    bool exiting;
    unsigned int id;

    // Output stream operator to print a PIDParticle
    friend std::ostream& operator<<(std::ostream& os, const PIDParticle& dt);
  };

  TSSelection(): _verbose(false) {}

  bool initialize(std::vector<std::string> input_files);

  bool run();

  bool analyze(gallery::Event* ev);

  bool finalize();

  void setVerbose(bool b) { _verbose = b; }

  void setOutputFile(TFile *f) { _fout = f; }

  // Set the producers for data products
  void setFluxWeightProducer(std::string s) { _fw_producer = s; }
  void setEventWeightProducer(std::string s) { _ew_producer = s; }
  void setMCTruthProducer(std::string s) { _mct_producer = s; }
  void setMCFluxProducer(std::string s) { _mcf_producer = s; }
  void setMCShowerProducer(std::string s) { _mcshw_producer = s; }
  void setMCTrackProducer(std::string s) { _mctrk_producer = s; }

  void setShowerEnergyResolution(float, bool);
  void setTrackEnergyResolution(float, bool);
  float nextShowerEnergyDistortion(float);
  float nextTrackEnergyDistortion(float);

  void setShowerAngleResolution(float, bool);
  void setTrackAngleResolution(float, bool);
  float nextShowerAngleDistortion(float);
  float nextTrackAngleDistortion(float);

  void setTrackShowerConfusion(float);
  bool nextTrackShowerConfusion();

  void setAcceptP(bool, int);
  void setAcceptNTrk(bool b) { _accept_ntrk = b; }

  void setNTrials(int n) { _n_trials = n; }

  // Set a numeric dataset ID, which is written into the tree as a tag
  void setDatasetID(int id) { _dataset_id = id; }

  // Utility function to test if a list of particles is 1lip
  // If EventType is not set, will accept any event enabled by selection.
  // Otherwise, will only accept events inside specified EventType.
  bool pass_selection(std::vector<PIDParticle>& p, int lpdg, EventType t=ANY); 

  // Apply track cuts
  static inline bool goodTrack(const sim::MCTrack& t, const simb::MCTruth& truth, float energy_distortion=0., float angle_distortion=0.) {
    return (!t.empty() &&
            tsutil::isFromNuVertex(truth, t) &&
            t.Process() == "primary" &&
            t.Start().E() - tsutil::get_pdg_mass(t.PdgCode()) + energy_distortion >= 60);
  }

  static inline bool goodTrack(bool isEmpty, bool isFromNuVertex, bool isPrimaryProcess, float energy) {
    return !isEmpty && isFromNuVertex && isPrimaryProcess && energy >= 60.;
  }

  // Apply shower cuts
  static inline bool goodShower(const sim::MCShower& s, const simb::MCTruth& truth, float energy_distortion=0., float angle_distortion=0.) {
    return (tsutil::isFromNuVertex(truth, s) &&
            s.Process() == "primary" &&
            (s.Start().E() - tsutil::get_pdg_mass(s.PdgCode())) + energy_distortion >= 30);
  }

  static inline bool goodShower(bool isFromNuVertex, bool isPrimaryProcess, float energy) {
    return isFromNuVertex && isPrimaryProcess && energy >= 30;
  }

  int get_nl(std::vector<PIDParticle>& p, int lpdg);
  int get_ntrk(std::vector<PIDParticle>& p);
  int get_np(std::vector<PIDParticle>& p);

  // A structure used to hold TTree output
  struct OutputData {
    OutputData() {
      weights = NULL;
    }
    int np;
    int n_trk;
    int nupdg;
    double enu;
    double q2;
    double w;
    double q0;
    double q3;
    int int_type;
    int int_mode;
    bool ccnc;
    double eccqe;
    std::vector<double> eps;
    std::vector<int> ppdgs;
    double elep;
    double thetalep;
    double philep;
    int lpdg;
    int lpid;
    double llen;
    bool lexit;
    double bnbweight;
    int dataset;
    std::map<std::string, std::vector<double> >* weights;
    unsigned event_number;
  };


  // structure to hold bokkeeping data
  struct HeaderData {
    // Data product producers
    std::string track_producer;
    std::string fw_producer;
    std::string ew_producer;
    std::string mct_producer;
    std::string mcf_producer;
    std::string mctrk_producer;
    std::string mcshw_producer;
    
    // Optionally set some energy resolution
    float shower_energy_resolution;
    bool shower_energy_by_percent;
    float track_energy_resolution;
    bool track_energy_by_percent;
    // and angle distortion
    float shower_angle_resolution;
    bool shower_angle_by_percent;
    float track_angle_resolution;
    bool track_angle_by_percent;
    
    // numbers of things
    int n_trials;
 
    // turn on/off different types of selections
    bool accept_0p;
    bool accept_1p;
    bool accept_np;
    bool accept_ntrk;
    
    // input files
    std::vector<std::string> input_files;  
    // event counts
    std::map<EventType, EventCounts> counts;
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
  std::map<EventType, EventCounts> _counts;

  // Optionally set some energy resolution
  float _shower_energy_resolution;
  bool _shower_energy_by_percent;
  std::normal_distribution<float> _shower_energy_distribution;
  float _track_energy_resolution;
  bool _track_energy_by_percent;
  std::normal_distribution<float> _track_energy_distribution;
  // and angular resolution
  float _shower_angle_resolution;
  bool _shower_angle_by_percent;
  std::normal_distribution<float> _shower_angle_distribution;
  float _track_angle_resolution;
  bool _track_angle_by_percent;
  std::normal_distribution<float> _track_angle_distribution;

  // track-shower confusion
  float _track_shower_confusion;
  std::bernoulli_distribution _track_shower_confusion_distribution;

  // random stuff
  std::mt19937 _gen;

  // number of times random stuff happens per event
  int _n_trials;

  // turn on/off different types of selections
  bool _accept_1p;
  bool _accept_0p;
  bool _accept_np;
  bool _accept_ntrk;

  // keep track of event index/number
  unsigned _event_number;

  bool _verbose;  //!< Print verbose output
  int _dataset_id;  //!< An arbitrary numeric ID
  OutputData* _data;  //!< Output data
  TTree* _tree;  //!< Output tree

  OutputData* _truth_data;  //!< Output data of truth
  TTree* _truth_data_tree;  //!< Output tree of truth

  TNtuple* _truthtree;
  TNtuple* _mectree;

  // input files
  std::vector<std::string> _input_files;

  // output file
  TFile* _fout;
  // whether to record truth level data
  bool _record_truth;
  bool _record_mec;

  TFile* _pdf_file;  //!< File containing dE/dx PDFs
  std::map<int, TH2F*> _trackdedxs;  // Track dE/dx distributions
};

}  // namespace galleryfmwk

#endif  // GALLERY_FMWK_TSSELECTION_H

