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
int main(int argv, char** argc) {
  if (argv < 3) {
    cout << "Must provide at least three input arguments: fout [fin]" << endl;
    return 1;
  }


  // All of the options
  int n_selections = 1;
  std::vector<int> dataset_id = {};

  std::vector<float> track_energy_distortion = {};
  bool track_energy_distortion_by_percent = false;

  std::vector<float> shower_energy_distortion = {};  
  bool shower_energy_distortion_by_percent = false;

  bool accept_np = false;
  bool accept_ntrack = false;

  cout << "Processing Runtime arguments" << endl;

  std::vector<TFile *> f_outs;
  //We have passed the input file as an argument to the function 
  vector<string> filename;
  for (int i = 1; i < argv; ++i) {
    if (strcmp(argc[i], "-o") == 0 || strcmp(argc[i], "--output_file") == 0) {
      i ++;
      while (i < argv && argc[i][0] != '-') {
        f_outs.push_back( new TFile(argc[i], "RECREATE") );
        cout << "Output File " << f_outs.size() << ": " << argc[i] << endl;
        i ++;
      }
    }
    if (strcmp(argc[i], "-n") == 0 || strcmp(argc[i], "--Nselections") == 0) {
      i ++;
      n_selections = stoi(argc[i]);
      cout << "N Selections: " << n_selections << endl;
      continue;
    }
    if (strcmp(argc[i], "-d") == 0 || strcmp(argc[i], "--datasetId") == 0) {
      while ( i+1 < argv && isNumber(argc[i+1]) ) {
        i ++;
        dataset_id.push_back( stoi(argc[i]) );
      }
      cout << "Dataset Id: " ;
      for (auto d: dataset_id) {
          cout << d << " ";
      }
      cout << endl;
      continue;
    }
    if (strcmp(argc[i], "--T_energy_distortion") == 0) {
      while( i+1 < argv && isNumber(argc[i+1]) ) {
        i ++;
        track_energy_distortion.push_back( stof(argc[i]) );
      }
      cout << "Track Energy Distortion: ";
      for (auto dist: track_energy_distortion) {
        cout << dist << " ";
      }
      cout << endl;
      continue;
    }
    if (strcmp(argc[i], "--T_edist_by_percent") == 0) {
      track_energy_distortion_by_percent = true;
      cout << "Track Energy Distortion By Percent"  << endl;
      continue;
    }
    if (strcmp(argc[i], "--S_energy_distortion") == 0) {
      while ( i+1 < argv && isNumber(argc[i+1]) ) {
        i ++;
        shower_energy_distortion.push_back( stof(argc[i]) );
      }
      cout << "Shower Energy Distortion: ";
      for (auto dist: shower_energy_distortion) {
        cout << dist << " ";
      }
      cout << endl;
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
    cout << "Input file: " << argc[i] << endl;
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

