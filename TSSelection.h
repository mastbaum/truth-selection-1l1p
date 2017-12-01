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

    // Output stream operator to print a PIDParticle
    friend std::ostream& operator<<(std::ostream& os, const PIDParticle& dt);
  };

  TSSelection(): _verbose(false) {}

  bool initialize();

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

  void setAcceptP(bool, int);
  void setAcceptNTrk(bool b) { _accept_ntrk = b; }

  // Set a numeric dataset ID, which is written into the tree as a tag
  void setDatasetID(int id) { _dataset_id = id; }

  // Utility function to test if a list of particles is 1lip
  bool is1lip(std::vector<PIDParticle>& p, int lpdg);

  // Apply track cuts
  static inline bool goodTrack(const sim::MCTrack& t, const simb::MCTruth& truth, float energy_distortion=0.) {
    return (!t.empty() &&
            tsutil::isFromNuVertex(truth, t) &&
            t.Process() == "primary" &&
            t.Start().E() - tsutil::get_pdg_mass(t.PdgCode()) + energy_distortion >= 60);
  }

  // Apply shower cuts
  static inline bool goodShower(const sim::MCShower& s, const simb::MCTruth& truth, float energy_distortion=0.) {
    return (tsutil::isFromNuVertex(truth, s) &&
            s.Process() == "primary" &&
            (s.Start().E() - tsutil::get_pdg_mass(s.PdgCode())) + energy_distortion >= 30);
  }

  // A structure used to hold TTree output
  struct OutputData {
    OutputData() {
      weights = NULL;
    }
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

  // Optionally set some energy resolution
  float _shower_energy_resolution;
  bool _shower_energy_by_percent;
  std::normal_distribution<float> _shower_energy_distribution;
  float _track_energy_resolution;
  bool _track_energy_by_percent;
  std::normal_distribution<float> _track_energy_distribution;
  // random stuff
  std::mt19937 _gen;

  // turn on/off different types of selections
  bool _accept_1p;
  bool _accept_np;
  bool _accept_ntrk;

  bool _verbose;  //!< Print verbose output
  int _dataset_id;  //!< An arbitrary numeric ID
  OutputData* _data;  //!< Output data
  TTree* _tree;  //!< Output tree

  TNtuple* _truthtree;
  TNtuple* _mectree;

  // output file
  TFile* _fout;

  TFile* _pdf_file;  //!< File containing dE/dx PDFs
  std::map<int, TH2F*> _trackdedxs;  // Track dE/dx distributions
};

}  // namespace galleryfmwk

#endif  // GALLERY_FMWK_TSSELECTION_H

