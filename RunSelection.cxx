#include <list>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TTree.h>

#include "TSUtil.h"
#include "TSConfig.h"
#include "TSSelection.h"

using namespace std;
using namespace art;

// copied from stack overflow
inline bool isInteger(const std::string & s)
{
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;

   char *p ;
   strtol(s.c_str(), &p, 10) ;

   return (*p == 0) ;
}

bool isNumber(const std::string& s)
{
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;

    char *end = 0;
    strtod(s.c_str(), &end);
    return end != s.c_str();
}

// stand alone executable for TSSelection
int main(int argc, char** argv) {
  if (argc < 3) {
    cout << "Must provide at least two input arguments: fout [fin]" << endl;
    return 1;
  }


  // All of the options
  int n_selections = 1;
  int n_trials = 1;
  std::vector<int> dataset_id = {};

  std::vector<float> track_energy_distortion = {};
  bool track_energy_distortion_by_percent = false;

  std::vector<float> shower_energy_distortion = {};  
  bool shower_energy_distortion_by_percent = false;

  bool drop_np = false;
  bool drop_ntrack = false;
  char *root_config_file_name = NULL;

  cout << "Processing Runtime arguments" << endl;

  std::vector<TFile *> f_outs;
  //We have passed the input file as an argument to the function 
  //
  // Argument parsing is fragile
  // if you input things incorrectly, it might not fail in a good way
  // probably use submitJobs.py to handle argument parsing
  //
  // If you set more than one selection, you should either pass in selection parameters 
  // like the energy resolution or once for every selection. 
  // I.e. the code fails unless energy_resolution.size() == 1 or energy_resolution.size() == n_selections
  //
  // takes in a list of input files as the arguments before any flags
  //
  // Since it is hard to parse particle id rates, things are set up to load from a config file. This config file
  // can be generated from gen_config.C.
  //
  // not all features in the selection are implemented as argument parsing
  vector<string> filename;
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output_file") == 0) {
      while (i+1 < argc && argv[i+1][0] != '-') {
        i ++;
        f_outs.push_back( new TFile(argv[i], "RECREATE") );
        cout << "Output File " << f_outs.size() << ": " << argv[i] << endl;
      }
      continue;
    }
    if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--Nselections") == 0) {
      i ++;
      n_selections = stoi(argv[i]);
      cout << "N Selections: " << n_selections << endl;
      continue;
    }
    if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--Ntrials") == 0) {
      i ++;
      n_trials = stoi(argv[i]);
      cout << "N Trials: " << n_trials << endl;
      continue;
    }
    if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--datasetId") == 0) {
      while ( i+1 < argc && isNumber(argv[i+1]) ) {
        i ++;
        dataset_id.push_back( stoi(argv[i]) );
      }
      cout << "Dataset Id: " ;
      for (auto d: dataset_id) {
          cout << d << " ";
      }
      cout << endl;
      continue;
    }
    if (strcmp(argv[i], "--T_energy_distortion") == 0) {
      while( i+1 < argc && isNumber(argv[i+1]) ) {
        i ++;
        track_energy_distortion.push_back( stof(argv[i]) );
      }
      cout << "Track Energy Distortion: ";
      for (auto dist: track_energy_distortion) {
        cout << dist << " ";
      }
      cout << endl;
      continue;
    }
    if (strcmp(argv[i], "--T_edist_by_percent") == 0) {
      track_energy_distortion_by_percent = true;
      cout << "Track Energy Distortion By Percent"  << endl;
      continue;
    }
    if (strcmp(argv[i], "--S_energy_distortion") == 0) {
      while ( i+1 < argc && isNumber(argv[i+1]) ) {
        i ++;
        shower_energy_distortion.push_back( stof(argv[i]) );
      }
      cout << "Shower Energy Distortion: ";
      for (auto dist: shower_energy_distortion) {
        cout << dist << " ";
      }
      cout << endl;
      continue;
    }
    if (strcmp(argv[i], "--S_edist_by_percent") == 0) {
      shower_energy_distortion_by_percent = true;
      cout << "Shower Energy Distortion By Percent"  << endl;
      continue;
    }
    if (strcmp(argv[i], "--drop_np") == 0) {
      drop_np = true;
      cout << "Drop Np" << endl;
      continue;
    }
    if (strcmp(argv[i], "--drop_ntrk") == 0) {
      drop_ntrack = true;
      cout << "Drop Ntrck" << endl;
      continue;
    }
    if (strcmp(argv[i], "--root_config_file") == 0) {
      i++;
      root_config_file_name = argv[i];
      continue;
    }
    filename.push_back(string(argv[i]));
    cout << "Input file: " << argv[i] << endl;
  }

  tsconfig::ConfigInfo *config = NULL;
  if (root_config_file_name != NULL) {
    config = tsconfig::ConfigInfo::load(root_config_file_name);
  }

  std::vector<galleryfmwk::TSSelection> selections(n_selections);

  cout << "Initialize" << endl;
  for (int i = 0; i < n_selections; i ++) {
    selections[i].setOutputFile(f_outs[i]); 
    std::vector<std::string> this_filename(filename);
    selections[i].initialize(this_filename);
  }

  for (int i = 0; i < n_selections; i ++) {
    selections[i].setFluxWeightProducer("eventweight");
    selections[i].setEventWeightProducer("mcweight");
    selections[i].setMCTruthProducer("generator");
    selections[i].setMCTrackProducer("mcreco");
    selections[i].setMCShowerProducer("mcreco");
    selections[i].setVerbose(true);
    selections[i].setNTrials(n_trials);
    selections[i].setAcceptP(!drop_np, 2);
    selections[i].setAcceptNTrk(!drop_ntrack);

    if (config != NULL) {
      assert(config->energy_range.size() == config->confusion_true_pdg.size());

      selections[i].setParticleIDEnergyRange(config->energy_range);
      for (unsigned j = 0; j < config->energy_range.size(); j++) {
        for (unsigned k = 0; k < config->confusion_true_pdg[j].size(); k++) {
          int true_pdgid = config->confusion_true_pdg[j][k];
          int test_pdgid = config->confusion_test_pdg[j][k];
          float id_rate = (float) config->confusion_id_rate[j][k];
          selections[i].addParticleIDRate(true_pdgid, test_pdgid, id_rate, (float) config->energy_range[j] - 1.);
          cout << "At ENERGY: " << config->energy_range[j] << " TRUE ID " << true_pdgid << " TEST ID " << test_pdgid << " RATE " << id_rate << endl;
        }
      }
      selections[i].checkParticleIDRates();
    }
   
    if (dataset_id.size() > 0) {
        selections[i].setDatasetID( dataset_id[i] );
    } 
    if (track_energy_distortion.size() > 0) {
      selections[i].setTrackEnergyResolution(track_energy_distortion[i], track_energy_distortion_by_percent);
    } 
    if (shower_energy_distortion.size() > 0) {
      selections[i].setShowerEnergyResolution(shower_energy_distortion[i], shower_energy_distortion_by_percent);
    }
  }

  cout << "Analyze" << endl;
  for (gallery::Event ev(filename) ; !ev.atEnd(); ev.next()) { 
    for (galleryfmwk::TSSelection &selection: selections) {
      selection.analyze(&ev);
    }
  }
  
  cout << "Finalize" << endl;
  for (galleryfmwk::TSSelection &selection: selections) {
    selection.finalize();
  }
}

