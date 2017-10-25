/**
 * Plot channel composition of CCQE energy distributions.
 *
 * It is assumed that (a) nue CC and (b) everything else come from different MC
 * with different POT scalings (hardcoded below, beware!).
 *
 * Usage:
 *
 *   $ root -l -x 'plot_eccqe.C("<input.root>")'
 *
 * (or compiled in ROOT with .L)
 *
 * Note: The POT scaling for the Monte Carlo samples (CCnue and !CCnue) are
 * hard-coded below and must be adjusted.
 *
 * A. Mastbaum <mastbaum@uchicago.edu>, 2017/10/25
 */

TH1F* get_channel(TTree* t, TString name, TString cut) {
  TH1F* h = new TH1F(name, ";E_{#nu}^{CCQE} (MeV);Events/100 MeV/6.6e20 POT", 25, 0, 2500);

  char sel[150];
  snprintf(sel, 150, "eccqe>>%s", name.Data());
  t->Draw(sel, cut, "goff");

  return h;  
}

void plot_eccqe(TString filename) {
  // POT for MC samples
  const float pot_ccnue = 8.47443e+21;
  const float pot_inclusive = 9.5184546e+20;

  gStyle->SetOptStat(1110);

  TFile* f = new TFile(filename);
  assert(f && f.IsOpen());
  TTree* data = (TTree*) f->Get("data");
  assert(data);

  // Histogram object naming:
  //
  //   h_T_QQ_CC_S
  //
  //  T  = e(lectron) for 1e1p, m(uon) for 1mu1p
  //  QQ = cc or nc
  //  CC = 00: QE, 01: RES, 02: DIS, 10: MEC
  //  S  = e for CCnue or i for BNB inclusive (sample)

  TH1F* h_m_cc_00_e = get_channel(data, "h_m_cc_00_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==0&&mode==0 )");
  TH1F* h_m_cc_01_e = get_channel(data, "h_m_cc_01_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==0&&mode==1 )");
  TH1F* h_m_cc_02_e = get_channel(data, "h_m_cc_02_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==0&&mode==2 )");
  TH1F* h_m_cc_10_e = get_channel(data, "h_m_cc_10_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==0&&mode==10)");

  TH1F* h_m_nc_00_e = get_channel(data, "h_m_nc_00_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==1&&mode==0 )");
  TH1F* h_m_nc_01_e = get_channel(data, "h_m_nc_01_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==1&&mode==1 )");
  TH1F* h_m_nc_02_e = get_channel(data, "h_m_nc_02_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==1&&mode==2 )");
  TH1F* h_m_nc_10_e = get_channel(data, "h_m_nc_10_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==1&&mode==10)");

  TH1F* h_m_cc_00_i = get_channel(data, "h_m_cc_00_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==0&&mode==0 )");
  TH1F* h_m_cc_01_i = get_channel(data, "h_m_cc_01_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==0&&mode==1 )");
  TH1F* h_m_cc_02_i = get_channel(data, "h_m_cc_02_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==0&&mode==2 )");
  TH1F* h_m_cc_10_i = get_channel(data, "h_m_cc_10_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==0&&mode==10)");

  TH1F* h_m_nc_00_i = get_channel(data, "h_m_nc_00_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==1&&mode==0 )");
  TH1F* h_m_nc_01_i = get_channel(data, "h_m_nc_01_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==1&&mode==1 )");
  TH1F* h_m_nc_02_i = get_channel(data, "h_m_nc_02_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==1&&mode==2 )");
  TH1F* h_m_nc_10_i = get_channel(data, "h_m_nc_10_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==13&&ccnc==1&&mode==10)");

  TH1F* h_e_cc_00_e = get_channel(data, "h_e_cc_00_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==0&&mode==0 )");
  TH1F* h_e_cc_01_e = get_channel(data, "h_e_cc_01_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==0&&mode==1 )");
  TH1F* h_e_cc_02_e = get_channel(data, "h_e_cc_02_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==0&&mode==2 )");
  TH1F* h_e_cc_10_e = get_channel(data, "h_e_cc_10_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==0&&mode==10)");

  TH1F* h_e_nc_00_e = get_channel(data, "h_e_nc_00_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==1&&mode==0 )");
  TH1F* h_e_nc_01_e = get_channel(data, "h_e_nc_01_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==1&&mode==1 )");
  TH1F* h_e_nc_02_e = get_channel(data, "h_e_nc_02_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==1&&mode==2 )");
  TH1F* h_e_nc_10_e = get_channel(data, "h_e_nc_10_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==1&&mode==10)");

  TH1F* h_e_cc_00_i = get_channel(data, "h_e_cc_00_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==0&&mode==0 )");
  TH1F* h_e_cc_01_i = get_channel(data, "h_e_cc_01_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==0&&mode==1 )");
  TH1F* h_e_cc_02_i = get_channel(data, "h_e_cc_02_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==0&&mode==2 )");
  TH1F* h_e_cc_10_i = get_channel(data, "h_e_cc_10_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==0&&mode==10)");

  TH1F* h_e_nc_00_i = get_channel(data, "h_e_nc_00_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==1&&mode==0 )");
  TH1F* h_e_nc_01_i = get_channel(data, "h_e_nc_01_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==1&&mode==1 )");
  TH1F* h_e_nc_02_i = get_channel(data, "h_e_nc_02_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==1&&mode==2 )");
  TH1F* h_e_nc_10_i = get_channel(data, "h_e_nc_10_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==11&&ccnc==1&&mode==10)");

  // Scale
  h_m_cc_00_i->Scale(6.6e20 / pot_inclusive);
  h_m_cc_01_i->Scale(6.6e20 / pot_inclusive);
  h_m_cc_02_i->Scale(6.6e20 / pot_inclusive);
  h_m_cc_10_i->Scale(6.6e20 / pot_inclusive);
  h_m_nc_00_i->Scale(6.6e20 / pot_inclusive);
  h_m_nc_01_i->Scale(6.6e20 / pot_inclusive);
  h_m_nc_02_i->Scale(6.6e20 / pot_inclusive);
  h_m_nc_10_i->Scale(6.6e20 / pot_inclusive);
  h_e_cc_00_i->Scale(6.6e20 / pot_inclusive);
  h_e_cc_01_i->Scale(6.6e20 / pot_inclusive);
  h_e_cc_02_i->Scale(6.6e20 / pot_inclusive);
  h_e_cc_10_i->Scale(6.6e20 / pot_inclusive);
  h_e_nc_00_i->Scale(6.6e20 / pot_inclusive);
  h_e_nc_01_i->Scale(6.6e20 / pot_inclusive);
  h_e_nc_02_i->Scale(6.6e20 / pot_inclusive);
  h_e_nc_10_i->Scale(6.6e20 / pot_inclusive);

  h_m_cc_00_e->Scale(6.6e20 / pot_ccnue);
  h_m_cc_01_e->Scale(6.6e20 / pot_ccnue);
  h_m_cc_02_e->Scale(6.6e20 / pot_ccnue);
  h_m_cc_10_e->Scale(6.6e20 / pot_ccnue);
  h_m_nc_00_e->Scale(6.6e20 / pot_ccnue);
  h_m_nc_01_e->Scale(6.6e20 / pot_ccnue);
  h_m_nc_02_e->Scale(6.6e20 / pot_ccnue);
  h_m_nc_10_e->Scale(6.6e20 / pot_ccnue);
  h_e_cc_00_e->Scale(6.6e20 / pot_ccnue);
  h_e_cc_01_e->Scale(6.6e20 / pot_ccnue);
  h_e_cc_02_e->Scale(6.6e20 / pot_ccnue);
  h_e_cc_10_e->Scale(6.6e20 / pot_ccnue);
  h_e_nc_00_e->Scale(6.6e20 / pot_ccnue);
  h_e_nc_01_e->Scale(6.6e20 / pot_ccnue);
  h_e_nc_02_e->Scale(6.6e20 / pot_ccnue);
  h_e_nc_10_e->Scale(6.6e20 / pot_ccnue);

  // Merge
  h_m_cc_00_i->Add(h_m_cc_00_e);
  h_m_cc_01_i->Add(h_m_cc_01_e);
  h_m_cc_02_i->Add(h_m_cc_02_e);
  h_m_cc_10_i->Add(h_m_cc_10_e);
  h_m_nc_00_i->Add(h_m_nc_00_e);
  h_m_nc_01_i->Add(h_m_nc_01_e);
  h_m_nc_02_i->Add(h_m_nc_02_e);
  h_m_nc_10_i->Add(h_m_nc_10_e);
  h_e_cc_00_i->Add(h_e_cc_00_e);
  h_e_cc_01_i->Add(h_e_cc_01_e);
  h_e_cc_02_i->Add(h_e_cc_02_e);
  h_e_cc_10_i->Add(h_e_cc_10_e);
  h_e_nc_00_i->Add(h_e_nc_00_e);
  h_e_nc_01_i->Add(h_e_nc_01_e);
  h_e_nc_02_i->Add(h_e_nc_02_e);
  h_e_nc_10_i->Add(h_e_nc_10_e);

  // Format
  h_m_cc_00_i->SetFillColor(kRed+1);
  h_m_cc_01_i->SetFillColor(kBlue+1);
  h_m_cc_02_i->SetFillColor(kGreen+1);
  h_m_cc_10_i->SetFillColor(kGray+1);
  h_m_nc_00_i->SetFillColor(kRed-7);
  h_m_nc_01_i->SetFillColor(kBlue-7);
  h_m_nc_02_i->SetFillColor(kGreen-7);
  h_m_nc_10_i->SetFillColor(kGray);
  h_e_cc_00_i->SetFillColor(kRed+1);
  h_e_cc_01_i->SetFillColor(kBlue+1);
  h_e_cc_02_i->SetFillColor(kGreen+1);
  h_e_cc_10_i->SetFillColor(kGray+1);
  h_e_nc_00_i->SetFillColor(kRed-7);
  h_e_nc_01_i->SetFillColor(kBlue-7);
  h_e_nc_02_i->SetFillColor(kGreen-7);
  h_e_nc_10_i->SetFillColor(kGray);

  // Legend
  TLegend* legend = new TLegend(0.65, 0.45, 0.88, 0.88);
  legend->AddEntry(h_m_cc_10_i, "CC MEC");
  legend->AddEntry(h_m_cc_02_i, "CC DIS");
  legend->AddEntry(h_m_cc_01_i, "CC RES");
  legend->AddEntry(h_m_cc_00_i, "CC QE");
  legend->AddEntry(h_m_nc_10_i, "NC MEC");
  legend->AddEntry(h_m_nc_02_i, "NC DIS");
  legend->AddEntry(h_m_nc_01_i, "NC RES");
  legend->AddEntry(h_m_nc_00_i, "NC QE");

  // MC statistical errors
  TH1F* h_m_i = get_channel(data, "h_m_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==13)");
  TH1F* h_m_e = get_channel(data, "h_m_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==13)");
  h_m_i->Sumw2();
  h_m_e->Sumw2();
  h_m_i->Scale(6.6e20 / pot_inclusive);
  h_m_e->Scale(6.6e20 / pot_ccnue);
  h_m_i->Add(h_m_e);
  h_m_i->SetLineWidth(2);

  TH1F* h_e_i = get_channel(data, "h_e_i", "bnbweight*(!(abs(nupdg)==12&&ccnc==0) && lpid==11)");
  TH1F* h_e_e = get_channel(data, "h_e_e", "bnbweight*( (abs(nupdg)==12&&ccnc==0) && lpid==11)");
  h_e_i->Sumw2();
  h_e_e->Sumw2();
  h_e_i->Scale(6.6e20 / pot_inclusive);
  h_e_e->Scale(6.6e20 / pot_ccnue);
  h_e_i->Add(h_e_e);
  h_e_i->SetLineWidth(2);

  // Stack
  THStack* hs_m = new THStack("mstack", ";E_{#nu}^{CCQE} (MeV);Events/100 MeV/6.6e20 POT");
  hs_m->Add(h_m_nc_00_i);
  hs_m->Add(h_m_nc_01_i);
  hs_m->Add(h_m_nc_02_i);
  hs_m->Add(h_m_nc_10_i);
  hs_m->Add(h_m_cc_00_i);
  hs_m->Add(h_m_cc_01_i);
  hs_m->Add(h_m_cc_02_i);
  hs_m->Add(h_m_cc_10_i);

  THStack* hs_e = new THStack("estack", ";E_{#nu}^{CCQE} (MeV);Events/100 MeV/6.6e20 POT");
  hs_e->Add(h_e_nc_00_i);
  hs_e->Add(h_e_nc_01_i);
  hs_e->Add(h_e_nc_02_i);
  hs_e->Add(h_e_nc_10_i);
  hs_e->Add(h_e_cc_00_i);
  hs_e->Add(h_e_cc_01_i);
  hs_e->Add(h_e_cc_02_i);
  hs_e->Add(h_e_cc_10_i);

  // Draw
  TCanvas* c1 = new TCanvas();
  hs_m->Draw("hist");
  h_m_i->Draw("e1x0 same");
  legend->Draw();
  c1->Update();
  c1->SaveAs("eccqe_1m1p.C");
  c1->SaveAs("eccqe_1m1p.pdf");
  c1->SaveAs("eccqe_1m1p.png");

  TCanvas* c2 = new TCanvas();
  hs_e->Draw("hist");
  h_e_i->Draw("e1x0 same");
  legend->Draw();
  c2->Update();
  c2->SaveAs("eccqe_1e1p.C");
  c2->SaveAs("eccqe_1e1p.pdf");
  c2->SaveAs("eccqe_1e1p.png");
}

