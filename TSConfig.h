/**
 * \file TSConfig.h
 */

#ifndef GALLERY_FMWK_TSCONFIG_H
#define GALLERY_FMWK_TSCONFIG_H

#include <assert.h>
#include <list>

#include <TFile.h>
#include <TTree.h>

namespace tsconfig {

  struct ConfigInfo {
    public:
    std::vector<float> energy_range;
    std::vector<std::vector<int>> confusion_true_pdg;
    std::vector<std::vector<int>> confusion_test_pdg;
    std::vector<std::vector<double>> confusion_id_rate;

    void save(const char *fname) {
      for (auto const &vector: confusion_id_rate) {
        double accumulator = 0.;
        for (double p: vector) {
          accumulator += p;
        }
        assert(abs(std::fmod(accumulator, 1.)) < 1e-4);
      }

      TFile *file = new TFile(fname, "RECREATE");
      file->cd();
      TTree *config = new TTree("data", "");
      config->Branch("energy_range", &this->energy_range);
      config->Branch("confusion_true_pdg", &this->confusion_true_pdg);
      config->Branch("confusion_test_pdg", &this->confusion_test_pdg);
      config->Branch("confusion_id_rate", &this->confusion_id_rate);

      config->Fill();
      config->Write();
    }

    ConfigInfo() : energy_range(), confusion_true_pdg(), confusion_test_pdg(), confusion_id_rate() {}

    ConfigInfo(std::list<std::list<std::tuple<int, int, double>>> confusion_input, std::list<float> energies_input) { 
      std::vector<std::vector<int>> true_pdg_vec;
      std::vector<std::vector<int>> test_pdg_vec; 
      std::vector<std::vector<double>> id_rate_vec;

      for (auto const &vector: confusion_input) {
        std::vector<int> confusion_true_pdg;
        std::vector<int> confusion_test_pdg;
        std::vector<double> confusion_rate_id;
        for (auto const &tuple: vector) {
          int pdg_true;
          int pdg_test;
          double id_rate;
          std::tie(pdg_true, pdg_test, id_rate) = tuple;
          confusion_true_pdg.push_back(pdg_true);
          confusion_test_pdg.push_back(pdg_test);
          confusion_rate_id.push_back(id_rate);
        }
        true_pdg_vec.push_back(confusion_true_pdg);
        test_pdg_vec.push_back(confusion_test_pdg);
        id_rate_vec.push_back(confusion_rate_id);
      }
      std::vector<float> energies_ret{ std::begin(energies_input), std::end(energies_input) };

      energy_range = energies_ret;
      confusion_true_pdg = true_pdg_vec;
      confusion_test_pdg = test_pdg_vec;
      confusion_id_rate = id_rate_vec;
    }

    static ConfigInfo *load(const char *fname) {
      TFile file(fname);
      file.cd();
      TTree* tree = (TTree*) file.Get("data");
      assert(tree && tree->GetEntries() > 0);

      ConfigInfo *ret = new ConfigInfo;

      auto energy_range_ref = &ret->energy_range;
      auto true_pdg_ref = &ret->confusion_true_pdg;
      auto test_pdg_ref = &ret->confusion_test_pdg;
      auto id_rate_ref = &ret->confusion_id_rate;

      tree->SetBranchAddress("energy_range", &energy_range_ref);
      tree->SetBranchAddress("confusion_true_pdg", &true_pdg_ref);
      tree->SetBranchAddress("confusion_test_pdg", &test_pdg_ref);
      tree->SetBranchAddress("confusion_id_rate", &id_rate_ref);

      tree->GetEntry(0);
      assert(ret->energy_range.size() > 0);
      assert(ret->confusion_true_pdg.size() > 0);
      for (auto const &vector: ret->confusion_id_rate) {
        double accumulator = 0.;
        for (double p: vector) {
          accumulator += p;
        }
        assert(abs(std::fmod(accumulator, 1.)) < 1e-4);
      }
      return ret;
    } 
  };

}  // namespace tsutil

#endif  // GALLERY_FMWK_TSCONFIG_H

