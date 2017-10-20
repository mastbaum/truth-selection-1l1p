/**
 * PID confusion matrix
 *
 * Print the true PDG composition of selected e/mu/p.
 *
 * Usage:
 *
 *   $ root -l -x 'confusion.C("<input.root>")'
 *
 *   (or compiled in ROOT with .L)
 *
 * A. Mastbaum <mastbaum@uchicago.edu>, 2017/09
 */

void confusion(TString filename) {
  gROOT->SetBatch(true);

  TFile* f = TFile::Open(filename);
  TTree* data = (TTree*) f->Get("data");

  // PID electrons
  float d1 = 100.0*data->Draw("eccqe","lpid==11&&lpdg== 11")/data->Draw("eccqe","lpid==11");
  float d2 = 100.0*data->Draw("eccqe","lpid==11&&lpdg==-11")/data->Draw("eccqe","lpid==11");
  float d3 = 100.0*data->Draw("eccqe","lpid==11&&lpdg== 22")/data->Draw("eccqe","lpid==11");

  std::cout << "e-/e-   " << d1 << std::endl;
  std::cout << "e-/e+   " << d2 << std::endl;
  std::cout << "e-/g    " << d3 << std::endl;
  //std::cout << d1+d2+d3 << std::endl;

  // PID muons
  float d4 = 100.0*data->Draw("eccqe","lpid==13&&lpdg==  13")/data->Draw("eccqe","lpid==13");
  float d5 = 100.0*data->Draw("eccqe","lpid==13&&lpdg==2212")/data->Draw("eccqe","lpid==13");
  float d6 = 100.0*data->Draw("eccqe","lpid==13&&lpdg== -13")/data->Draw("eccqe","lpid==13");
  float d7 = 100.0*data->Draw("eccqe","lpid==13&&lpdg== 211")/data->Draw("eccqe","lpid==13");
  float d8 = 100.0*data->Draw("eccqe","lpid==13&&lpdg==-211")/data->Draw("eccqe","lpid==13");
  float dk = 100.0*data->Draw("eccqe","lpid==13&&lpdg== 321")/data->Draw("eccqe","lpid==13");

  std::cout << "mu-/mu- " << d4 << std::endl;
  std::cout << "mu-/p   " << d5 << std::endl;
  std::cout << "mu-/mu+ " << d6 << std::endl;
  std::cout << "mu-/pi+ " << d7 << std::endl;
  std::cout << "mu-/pi- " << d8 << std::endl;
  std::cout << "mu-/k+  " << dk << std::endl;
  //std::cout << d4+d5+d6+d7+d8+dk << std::endl;

  // PID protons
  float d9 = 100.0*data->Draw("eccqe","ppdg==  13")/data->Draw("eccqe","1.0");
  float da = 100.0*data->Draw("eccqe","ppdg==2212")/data->Draw("eccqe","1.0");
  float db = 100.0*data->Draw("eccqe","ppdg== -13")/data->Draw("eccqe","1.0");
  float dc = 100.0*data->Draw("eccqe","ppdg== 211")/data->Draw("eccqe","1.0");
  float dd = 100.0*data->Draw("eccqe","ppdg==-211")/data->Draw("eccqe","1.0");
  float de = 100.0*data->Draw("eccqe","ppdg== 321")/data->Draw("eccqe","1.0");
  float df = 100.0*data->Draw("eccqe","ppdg==2112")/data->Draw("eccqe","1.0");

  std::cout << "p/mu-   " << d9 << std::endl;
  std::cout << "p/p     " << da << std::endl;
  std::cout << "p/mu+   " << db << std::endl;
  std::cout << "p/pi+   " << dc << std::endl;
  std::cout << "p/pi-   " << dd << std::endl;
  std::cout << "p/k+    " << de << std::endl;
  std::cout << "p/n     " << df << std::endl;
  //std::cout << d9+da+db+dc+dd+de+df << std::endl;
}

