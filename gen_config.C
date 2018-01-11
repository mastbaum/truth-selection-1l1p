#include "TSConfig.h"

void gen_config() {
  std::list<float> energy_range = {0., 100., 200.};
  std::list<std::list<std::tuple<int, int, double>>> confusion_entries = 
    {
        {std::tuple<int, int, double>(0, 0, 0.5), std::tuple<int, int, double>(0, 1, 0.5), std::tuple<int, int, double>(1, 1, 1.)} ,
        {std::tuple<int, int, double>(0, 0, 0.5), std::tuple<int, int, double>(0, 1, 0.5), std::tuple<int, int, double>(1, 1, 1.)} ,
    };
  tsconfig::ConfigInfo config = tsconfig::ConfigInfo(confusion_entries, energy_range);

  config.save("config.root");
}

