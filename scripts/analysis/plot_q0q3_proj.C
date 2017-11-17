/**
 * Plot q0 vs. q3
 *
 * It is assumed that (a) nue CC and (b) everything else come from different MC
 * with different POT scalings (hardcoded below, beware!).
 *
 * Usage:
 *
 *   $ root -l -x 'plot_q0q3.C("<input.root>")'
 *
 * (or compiled in ROOT with .L)
 *
 * Note: The POT scaling for the Monte Carlo samples (CCnue and !CCnue) are
 * hard-coded below and must be adjusted.
 *
 * A. Mastbaum <mastbaum@uchicago.edu>, 2017/10
 */

void plot_q0q3_proj(TString filename) {
  const float pot_ccnue = 8.47443e+21;
  const float pot_inclusive = 9.5184546e+20;

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas();

  TFile* f = new TFile(filename);
  assert(f && f.IsOpen());
  TTree* data = (TTree*) f->Get("data");
  assert(data);

  data->SetLineColor(kWhite);
  data->SetMarkerStyle(6);

  data->SetMarkerColor(kBlue-7);
  data->Draw("q0-(sqrt(0.938**2+q3**2)-0.938):q3>>h_m_i(100,0,1.2,100,-0.6,0.6)",
             "bnbweight*(lpid==13&&!(abs(nupdg)==12&&ccnc==0&&mode==0))");
  data->Draw("q0-(sqrt(0.938**2+q3**2)-0.938):q3>>h_m_e(100,0,1.2,100,-0.6,0.6)",
             "bnbweight*(lpid==13&& (abs(nupdg)==12&&ccnc==0&&mode==0))");

  data->SetMarkerColor(kRed);
  data->Draw("q0-(sqrt(0.938**2+q3**2)-0.938):q3>>h_e_i(100,0,1.2,100,-0.6,0.6)",
             "bnbweight*(lpid==11&&!(abs(nupdg)==12&&ccnc==0&&mode==0))");
  data->Draw("q0-(sqrt(0.938**2+q3**2)-0.938):q3>>h_e_e(100,0,1.2,100,-0.6,0.6)",
             "bnbweight*(lpid==11&& (abs(nupdg)==12&&ccnc==0&&mode==0))");

  TH2F* h_m_i = (TH2F*) gDirectory->Get("h_m_i");
  TH2F* h_m_e = (TH2F*) gDirectory->Get("h_m_e");
  TH2F* h_e_i = (TH2F*) gDirectory->Get("h_e_i");
  TH2F* h_e_e = (TH2F*) gDirectory->Get("h_e_e");

  h_m_i->Scale(6.6e20 / pot_inclusive);
  h_m_e->Scale(6.6e20 / pot_ccnue);
  h_e_i->Scale(6.6e20 / pot_inclusive);
  h_e_e->Scale(6.6e20 / pot_ccnue);

  h_m_i->Add(h_m_e);
  h_e_i->Add(h_e_e);

  h_e_i->Scale(5000.0 / h_e_i->Integral());
  h_m_i->Scale(5000.0 / h_m_i->Integral());

  h_e_i->SetMarkerColor(kRed);
  h_e_i->SetMarkerStyle(7);

  h_m_i->SetMarkerColor(kBlue-7);
  h_m_i->SetMarkerStyle(7);

  h_e_i->SetTitle(";True three-momentum transfer (GeV);True energy transfer (GeV)");
  h_e_i->Draw();
  h_m_i->Draw("same");

  //TF1* f1 = new TF1("f_w938","sqrt(0.938**2 + x**2) - 0.938" , 0.0, 1.2);
  //f1->SetLineColor(kBlack);
  //f1->SetLineWidth(1);
  //f1->Draw("same");

  TLegend* l = new TLegend(0.14, 0.65, 0.75, 0.87);
  //l->AddEntry(q2lines[0], "Q^{2} = 0.2 to 1.0 GeV^{2}");
  l->AddEntry(h_m_i, "1#mu1p selection");
  l->AddEntry(h_e_i, "1e1p selection");
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->Draw();

  c1->Update();
  //c1->SaveAs("q0q3.C");
  //c1->SaveAs("q0q3.pdf");
  //c1->SaveAs("q0q3.png");
}

