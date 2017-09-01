void plot() {
  TFile* f = TFile::Open("sel_inc_short.root");
  TNtuple* data = (TNtuple*) f->Get("data");

  THStack* hs = new THStack("hs", "");
  data->Draw("eccqe>>h0cc( 25,0,2500)", "lpid==11&&ccnc==0&&mode==0", "goff");
  data->Draw("eccqe>>h1cc( 25,0,2500)", "lpid==11&&ccnc==0&&mode==1", "goff");
  data->Draw("eccqe>>h2cc( 25,0,2500)", "lpid==11&&ccnc==0&&mode==2", "goff");
  data->Draw("eccqe>>h10cc(25,0,2500)", "lpid==11&&ccnc==0&&mode==10", "goff");

  //data->Draw("eccqe>>h0nc( 25,0,2500)", "lpid==13&&mode== 0&&ccnc==1", "goff");
  //data->Draw("eccqe>>h1nc( 25,0,2500)", "lpid==13&&mode== 1&&ccnc==1", "goff");
  //data->Draw("eccqe>>h2nc( 25,0,2500)", "lpid==13&&mode== 2&&ccnc==1", "goff");
  //data->Draw("eccqe>>h10nc(25,0,2500)", "lpid==13&&mode==10&&ccnc==1", "goff");

  TH1F* h0cc  = (TH1F*) gDirectory->Get("h0cc");
  TH1F* h1cc  = (TH1F*) gDirectory->Get("h1cc");
  TH1F* h2cc  = (TH1F*) gDirectory->Get("h2cc");
  TH1F* h10cc = (TH1F*) gDirectory->Get("h10cc");

  //TH1F* h0nc  = (TH1F*) gDirectory->Get("h0nc");
  //TH1F* h1nc  = (TH1F*) gDirectory->Get("h1nc");
  //TH1F* h2nc  = (TH1F*) gDirectory->Get("h2nc");
  //TH1F* h10nc = (TH1F*) gDirectory->Get("h10nc");

  h0cc->SetFillColor(kRed+1);
  h1cc->SetFillColor(kBlue+1);
  h2cc->SetFillColor(kGreen+1);
  h10cc->SetFillColor(kGray+1);

  //h0nc->SetFillColor(kRed+2);
  //h1nc->SetFillColor(kBlue+2);
  //h2nc->SetFillColor(kGreen+2);
  //h10nc->SetFillColor(kGray+2);

  hs->Add(h0cc);
  hs->Add(h10cc);
  hs->Add(h1cc);
  hs->Add(h2cc);

  //hs->Add(h0nc);
  //hs->Add(h10nc);
  //hs->Add(h1nc);
  //hs->Add(h2nc);

  TLegend* l = new TLegend(0.6, 0.6, 0.88, 0.88);
  l->AddEntry(h0cc,  "QE");
  l->AddEntry(h1cc,  "RES");
  l->AddEntry(h2cc,  "DIS");
  l->AddEntry(h10cc, "MEC");

  //l->AddEntry(h0nc, "NCQE");
  //l->AddEntry(h1nc, "NCRES");
  //l->AddEntry(h2nc, "NCDIS");
  //l->AddEntry(h10nc, "NCMEC");

  //hs->SetTitle("1#mu1p selection;E_{#nu}^{CCQE} (MeV);Events");
  hs->SetTitle(";E_{#nu}^{CCQE} (MeV);Events");
  hs->Draw();
  l->Draw();

  gPad->Update();
  //TPaveText* title = (TPaveText*) gPad->GetPrimitive("title");
  //title->SetBorderSize(0);
  //title->SetFillColor(kWhite);
  //title->SetTextFont(132);
  gPad->Update();
  gPad->Update();
}

