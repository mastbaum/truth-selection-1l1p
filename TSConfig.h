/**
 * \file TSConfig.h
 */

#ifndef GALLERY_FMWK_TSCONFIG_H
#define GALLERY_FMWK_TSCONFIG_H

#include <assert.h>
#include <list>

namespace tsconfig {

  struct ConfigInfo {
    public:
    std::list<float> energy_range;
    std::list<std::list<std::tuple<int, int, float>>> confusion_entries;

    void save(const char *fname) {
      TFile file(fname, "RECREATE");
      file.cd();
      TTree *config = new TTree("config", "");
      TBranch *energy_range = config->Branch("energy_range", &energy_range);
      TBranch *confusion_entries = config->Branch("confusion_entries", &confusion_entries);
      config->Fill();
      config->Write();
    }

    std::vector<std::vector<std::tuple<int, int, float>>> confusionEntries() {
      std::vector<std::vector<std::tuple<int, int, float>>> ret;
      for (auto const &list: confusion_entries) {
         std::vector<std::tuple<int, int, float>> data{ std::begin(list), std::end(list) };
         ret.push_back(data);
      }
      return ret;
    }

    std::vector<float> energyRange() {
      std::vector<float> ret{ std::begin(energy_range), std::end(energy_range) };
      return ret;
    }

    static ConfigInfo *load(const char *fname) {
      TFile file(fname);
      file.cd();
      TTree* tree = (TTree*) file.Get("config");
      assert(tree && tree->GetEntries() > 0);
      ConfigInfo *ret = new ConfigInfo;
      tree->SetBranchAddress("energy_range", &ret->energy_range);
      tree->SetBranchAddress("confusion_entries", &ret->confusion_entries);
      tree->GetEntry(0);
      assert(ret->energy_range.size() > 0);
      assert(ret->confusion_entries.size() > 0);
      return ret;
    } 
  };

}  // namespace tsutil

#endif  // GALLERY_FMWK_TSCONFIG_H

