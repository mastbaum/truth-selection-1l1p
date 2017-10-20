/**
 * Add together nues from the nue-only sample, and numus from the
 * inclusive sample, to build up a high-stats combined sample.
 *
 * We take events with "abs(nupdg)==12&&ccnc==0" -- CC nue and nuebar
 * interactions -- from the CCnue file, and everything but that subset
 * from the inclusive BNB file. This exclusivity lets us easily do POT
 * scaling based on that event specification.
 *
 * Usage:
 *
 *   root -l -x -q add.C
 *
 * N.B. The input and output filenames are hard-coded at the top of the
 * script.
 *
 * A. Mastbaum <mastbaum@uchicago.edu>, 2017/09/08
 */

{
  // CCnue input
  TString f_ccnue = "../sel_ccnue_take2.root";

  // BNB inclusive input
  TString f_bnb = "../sel_inclusive_mcc8a_try2.root";

  // Output
  TString f_out = "sel_merged_mcc8.root";


  // Read in the input files
  TFile* fe = TFile::Open(f_ccnue);
  TTree* te = (TTree*) fe->Get("data");
  std::cout << "ccnue all: " << te->GetEntries() << std::endl;

  TFile* fi2 = TFile::Open(f_bnb);
  TTree* ti2 = (TTree*) fi2->Get("data");
  std::cout << "inc all: " << ti2->GetEntries() << std::endl;

  // Construct the output file
  TFile* fo = TFile::Open(f_out, "recreate");

  // From the CCnue sample, grab CC nue and nuebar
  TTree* tccnue = te->CopyTree("abs(nupdg)==12&&ccnc==0");
  std::cout << "ccnue: " << tccnue->GetEntries() << std::endl;

  // From the BNB inclusive sample, grab everything else
  TTree* tother2 = ti2->CopyTree("!(abs(nupdg)==12&&ccnc==0)");
  std::cout << "other: " << tother2->GetEntries() << std::endl;

  TList* list = new TList;
  list->Add(tccnue);
  list->Add(tother2);
  TTree* tc = TTree::MergeTrees(list);
  tc->SetName("data");

  std::cout << "output: "  << tc->GetEntries() << std::endl;

  tc->Write();
  fo->Close();
}

