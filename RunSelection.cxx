#include <experimental/optional>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TTree.h>

#include "TSUtil.h"
#include "TSSelection.h"

using namespace std;
using namespace art;

// stand alone executable for TSSelection
int main(int argv, char** argc) {
  if (argv < 3) {
    cout << "Must provide at least three input arguments: fout [fin]" << endl;
    return 1;
  }


  // All of the options
  std::experimental::optional<int> maybe_dataset_id = {};

  std::experimental::optional<float> maybe_track_energy_distortion = {};
  bool track_energy_distortion_by_percent = false;

  std::experimental::optional<float> maybe_shower_energy_distortion = {};  
  bool shower_energy_distortion_by_percent = false;

  bool accept_np = false;
  bool accept_ntrack = false;

  cout << "Processing Runtime arguments" << endl;

  TFile f_out(argc[1], "NEW");
  //We have passed the input file as an argument to the function 
  vector<string> filename;
  for (int i = 2; i < argv; ++i) {
    if (strcmp(argc[i], "-d") == 0 || strcmp(argc[i], "--datasetId") == 0) {
      i ++;
      maybe_dataset_id = stoi(argc[i]);
      cout << "Dataset Id: " << maybe_dataset_id.value() << endl;
      continue;
    }
    if (strcmp(argc[i], "--T_energy_distortion") == 0) {
      i ++;
      maybe_track_energy_distortion = stof(argc[i]);
      cout << "Track Energy Distortion: " << maybe_track_energy_distortion.value() << endl;
      continue;
    }
    if (strcmp(argc[i], "--T_edist_by_percent") == 0) {
      track_energy_distortion_by_percent = true;
      cout << "Track Energy Distortion By Percent"  << endl;
      continue;
    }
    if (strcmp(argc[i], "--S_energy_distortion") == 0) {
      i ++;
      maybe_shower_energy_distortion = stof(argc[i]);
      cout << "Shower Energy Distortion: " << maybe_shower_energy_distortion.value() << endl;
      continue;
    }
    if (strcmp(argc[i], "--S_edist_by_percent") == 0) {
      shower_energy_distortion_by_percent = true;
      cout << "Shower Energy Distortion By Percent"  << endl;
      continue;
    }
    if (strcmp(argc[i], "--accept_np") == 0) {
      accept_np = true;
      cout << "Accept Np" << endl;
      continue;
    }
    if (strcmp(argc[i], "--accept_ntrk") == 0) {
      accept_ntrack = true;
      cout << "Accept Ntrck" << endl;
      continue;
    }
    filename.push_back(string(argc[i]));
  }

  galleryfmwk::TSSelection selection;
  selection.setFluxWeightProducer("eventweight");
  selection.setEventWeightProducer("mcweight");
  selection.setMCTruthProducer("generator");
  selection.setMCTrackProducer("mcreco");
  selection.setMCShowerProducer("mcreco");
  selection.setVerbose(true);

  selection.setOutputFile(&f_out); 

  if (bool(maybe_dataset_id)) {
    selection.setDatasetID(maybe_dataset_id.value());
  }
  if (bool(maybe_track_energy_distortion)) {
    selection.setTrackEnergyResolution(maybe_track_energy_distortion.value(), track_energy_distortion_by_percent);
  } 
  if (bool(maybe_shower_energy_distortion)) {
    selection.setShowerEnergyResolution(maybe_shower_energy_distortion.value(), shower_energy_distortion_by_percent);
  }

  cout << "Initialized" << endl;
  selection.initialize();

  cout << "Analyze" << endl;
  for (gallery::Event ev(filename) ; !ev.atEnd(); ev.next()) {
    selection.analyze(&ev);
  }
  
  cout << "Finalize" << endl;
  selection.finalize();
}

