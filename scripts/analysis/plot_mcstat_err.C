/**
 * Plot errors for Monte Carlo statistics.
 *
 * It is assumed that (a) nue CC and (b) everything else come from different MC
 * with different POT scalings (hardcoded below, beware!).
 *
 * Usage:
 *
 *   $ root -l -x 'plot_mcstat_err.C("<input.root>")'
 *
 * (or compiled in ROOT with .L)
 *
 * Note: The POT scaling for the Monte Carlo samples (CCnue and !CCnue) are
 * hard-coded below and must be adjusted.
 *
 * A. Mastbaum <mastbaum@uchicago.edu>, 2017/10
 */

TH1F* get_channel(TTree* t, TString name, TString cut) {
  TH1F* h = new TH1F(name, ";E_{#nu}^{CCQE} (MeV);Events/100 MeV/6.6e20 POT", 25, 0, 2500);

  char sel[150];
  snprintf(sel, 150, "eccqe>>%s", name.Data());
  t->Draw(sel, cut, "goff");

  return h;  
}

void plot_mcstat_err(TString filename) {
  const float pot_ccnue = 8.47443e+21;
  const float pot_inclusive = 9.5184546e+20;

  gStyle->SetOptStat(1110);

  TFile* f = new TFile(filename);
  assert(f && f.IsOpen());
  TTree* data = (TTree*) f->Get("data");
  assert(data);

  // MC statistical errors
  TH1F* h_m_i = get_channel(data, "h_m_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==13)");
  TH1F* h_m_e = get_channel(data, "h_m_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==13)");
  h_m_i->Sumw2();
  h_m_e->Sumw2();
  h_m_i->Scale(6.6e20 / pot_inclusive);
  h_m_e->Scale(6.6e20 / pot_ccnue);
  h_m_i->Add(h_m_e);
  //h_m_i->SetLineWidth(2);

  TH1F* h_e_i = get_channel(data, "h_e_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==11)");
  TH1F* h_e_e = get_channel(data, "h_e_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==11)");
  h_e_i->Sumw2();
  h_e_e->Sumw2();
  h_e_i->Scale(6.6e20 / pot_inclusive);
  h_e_e->Scale(6.6e20 / pot_ccnue);
  h_e_i->Add(h_e_e);
  //h_e_i->SetLineWidth(2);


  //data->Draw("eccqe>>h1(25,0,2500)","lpid==11&&!(abs(nupdg)==12&&ccnc==0)");
  //data->Draw("eccqe>>h2(25,0,2500)","lpid==11&& (abs(nupdg)==12&&ccnc==0)");
  //TH1F* h1 = (TH1F*) gDirectory->Get("h1");
  //TH1F* h2 = (TH1F*) gDirectory->Get("h2");
  //h1->Sumw2();
  //h2->Sumw2();

  //// POT scaling
  //h2->Scale(6.6e20 / pot_inclusive);
  //h1->Scale(6.6e20 / pot_ccnue);
  ////h1->Scale(6.6e20 / (2.809038e+19 + 1.1386693e+20));

  TCanvas* c1 = new TCanvas();
  h_m_i->SetTitle(";E_{#nu}^{CCQE};Events/100 MeV/6.6e20 POT");
  h_m_i->Draw("e2");
  h_m_i->SetFillColor(632-7);
  h_m_i->SetMarkerColor(632-7);
  c1->Update();
  c1->SaveAs("mcstat_errors_1m1p.C");
  c1->SaveAs("mcstat_errors_1m1p.pdf");
  c1->SaveAs("mcstat_errors_1m1p.png");
  std::cout << "1m1p:" << std::endl;
  for (int i=0; i<h_m_i->GetNbinsX(); i++) {
    std::cout << i << " " << h_m_i->GetBinCenter(i) << " " << h_m_i->GetBinError(i) / h_m_i->GetBinContent(i) << std::endl;
  }

  TCanvas* c2 = new TCanvas();
  h_e_i->SetTitle(";E_{#nu}^{CCQE};Events/100 MeV/6.6e20 POT");
  h_e_i->Draw("e2");
  h_e_i->SetFillColor(632-7);
  h_e_i->SetMarkerColor(632-7);
  c2->Update();
  c2->SaveAs("mcstat_errors_1e1p.C");
  c2->SaveAs("mcstat_errors_1e1p.pdf");
  c2->SaveAs("mcstat_errors_1e1p.png");
  std::cout << "1e1p:" << std::endl;
  for (int i=0; i<h_e_i->GetNbinsX(); i++) {
    std::cout << i << " " << h_e_i->GetBinCenter(i) << " " << h_e_i->GetBinError(i) / h_e_i->GetBinContent(i) << std::endl;
  }
}

