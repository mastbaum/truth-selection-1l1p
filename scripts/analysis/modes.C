/**
 * Plot true interaction modes.
 *
 * Draws cos(theta_lep) vs. eccqe and eccqe for CC QE/RES/DIS/COH/MEC.
 *
 * Usage:
 *
 *   $ root -l -x 'modes.C("<input.root>")'
 *
 *   (or compiled in ROOT with .L)
 *
 * A. Mastbaum <mastbaum@uchicago.edu>, 2017/09
 */

void modes(TString filename) {
  TFile* f = TFile::Open(filename);
  TTree* d = (TTree*) f->Get("data");

  d->SetMarkerStyle(7);

  d->SetMarkerColor(kRed);
  d->Draw("cos(thetalep):eccqe>>h00", "ep<1e3&&ccnc==0&&mode==0");
  TH2F* h00 = (TH2F*) gDirectory->FindObject("h00");
  h00->Scale(1.0 / h00->Integral());
  h00->SetTitle(";E_{#nu}^{CCQE};cos#theta_{lep}");

  d->SetMarkerColor(kBlue);
  d->Draw("cos(thetalep):eccqe>>h01", "ep<1e3&&ccnc==0&&mode==1", "same");
  TH2F* h01 = (TH2F*) gDirectory->FindObject("h01");
  h01->Scale(1.0 / h01->Integral());

  d->SetMarkerColor(kGreen);
  d->Draw("cos(thetalep):eccqe>>h02", "ep<1e3&&ccnc==0&&mode==2", "same");
  TH2F* h02 = (TH2F*) gDirectory->FindObject("h02");
  h02->Scale(1.0 / h02->Integral());

  d->SetMarkerColor(kGray);
  d->Draw("cos(thetalep):eccqe>>h10", "ep<1e3&&ccnc==0&&mode==10", "same");
  TH2F* h10 = (TH2F*) gDirectory->FindObject("h10");
  h10->Scale(1.0 / h10->Integral());

  TLegend* l2 = new TLegend(0.6, 0.6, 0.85, 0.85);
  l2->AddEntry(h00, "CC QE");
  l2->AddEntry(h01, "CC RES");
  l2->AddEntry(h02, "CC COH");
  l2->AddEntry(h10, "CC MEC");
  l2->Draw();

  TCanvas* c2 = new TCanvas();
  gStyle->SetOptStat(0);
  d->SetMarkerColor(kBlack);
  d->SetMarkerStyle(0);

  d->SetLineColor(kRed);
  d->Draw("eccqe>>hc00(25,0,3000)", "lpid==13&&ccnc==0&&mode==0", "hist");
  TH1F* hc00 = (TH1F*) gDirectory->FindObject("hc00");
  hc00->Scale(1.0 / hc00->Integral());
  hc00->SetLineWidth(2);
  hc00->SetTitle(";E_{#nu}^{CCQE} (MeV);Probability/bin");
  hc00->GetYaxis()->SetRangeUser(0, 0.12);

  d->SetLineColor(kBlue);
  d->Draw("eccqe>>hc01(25,0,3000)", "lpid==13&&ccnc==0&&mode==1", "same hist");
  TH1F* hc01 = (TH1F*) gDirectory->FindObject("hc01");
  hc01->Scale(1.0 / hc01->Integral());
  hc01->SetLineWidth(2);

  d->SetLineColor(kGreen+1);
  d->Draw("eccqe>>hc02(25,0,3000)", "lpid==13&&ccnc==0&&mode==2", "same hist");
  TH1F* hc02 = (TH1F*) gDirectory->FindObject("hc02");
  hc02->Scale(1.0 / hc02->Integral());
  hc02->SetLineWidth(2);

  d->SetLineColor(kGray);
  d->Draw("eccqe>>hc10(25,0,3000)", "lpid==13&&ccnc==0&&mode==10", "same hist");
  TH1F* hc10 = (TH1F*) gDirectory->FindObject("hc10");
  hc10->Scale(1.0 / hc10->Integral());
  hc10->SetLineWidth(2);

  //d->SetLineColor(kOrange);
  //d->Draw("eccqe>>hco(25,0,3000)", "lpid==11&&mode!=10&&mode>2", "same hist");
  //TH1F* hco = (TH1F*) gDirectory->FindObject("hco");
  //hco->Scale(1.0 / hco->Integral());
  //hco->SetLineWidth(2);

  TLegend* l = new TLegend(0.6, 0.6, 0.85, 0.85);
  l->AddEntry(hc00, "CC QE");
  l->AddEntry(hc01, "CC RES");
  l->AddEntry(hc02, "CC COH");
  l->AddEntry(hc10, "CC MEC");
  //l->AddEntry(hco, "Other");
  l->Draw();
}

