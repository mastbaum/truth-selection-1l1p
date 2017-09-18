/**
 * Add together nues from the nue-only sample, and numus from the
 * inclusive sample, to build up a high-stats combined sample.
 *
 * A. Mastbaum <mastbaum@uchicago.edu>, 2017/09/08
 */

{
  TFile* fe = TFile::Open("../sel_e.root");
  TTree* te = (TTree*) fe->Get("data");
  std::cout << "nue all: " << te->GetEntries() << std::endl;

  TFile* fi = TFile::Open("../sel_inclusive.root");
  TTree* ti = (TTree*) fi->Get("data");
  std::cout << "inc all: " << ti->GetEntries() << std::endl;

  TFile* fo = TFile::Open("./sel_e_mu.root", "recreate");

  TTree* teel = te->CopyTree("lpid==11");
  std::cout << "nue e: "   << teel->GetEntries() << std::endl;
  TTree* timu = ti->CopyTree("lpid==13");
  std::cout << "inc mu: "  << timu->GetEntries() << std::endl;

  TList* list = new TList;
  list->Add(teel);
  list->Add(timu);
  TTree* tc = TTree::MergeTrees(list);
  tc->SetName("data");

  std::cout << "output: "  << tc->GetEntries() << std::endl;

  tc->Write();
  fo->Close();
}

