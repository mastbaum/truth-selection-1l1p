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

void plot_mcstat_err(TString filename) {
  gStyle->SetOptStat(1110);
  TCanvas* c1 = new TCanvas();

  TFile* f = new TFile(filename);
  assert(f && f.IsOpen());
  TTree* data = (TTree*) f->Get("data");
  assert(data);

  data->Draw("eccqe>>h1(25,0,2500)","lpid==11&&!(abs(nupdg)==12&&ccnc==0)");
  data->Draw("eccqe>>h2(25,0,2500)","lpid==11&& (abs(nupdg)==12&&ccnc==0)");
  TH1F* h1 = (TH1F*) gDirectory->Get("h1");
  TH1F* h2 = (TH1F*) gDirectory->Get("h2");
  h1->Sumw2();
  h2->Sumw2();

  // POT scaling
  h2->Scale(6.6e20 / 8.47443e+21);
  h1->Scale(6.6e20 / (2.809038e+19 + 1.1386693e+20));

  h2->Add(h1);
  h2->SetTitle(";E_{#nu}^{CCQE};Events/100 MeV/6.6e20 POT");
  h2->Draw("e2");
  h2->SetFillColor(632-7);
  h2->SetMarkerColor(632-7);
  c1->Update();
  c1->SaveAs("mcstat_errors_1e1p.C");
  c1->SaveAs("mcstat_errors_1e1p.pdf");
  c1->SaveAs("mcstat_errors_1e1p.png");
}

