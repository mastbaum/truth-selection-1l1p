/**
 * Plot bin-to-bin correlations.
 *
 * Produce a matrix-like plot where each cell is a covariance matrix element,
 * but shows the actual distribution of bin contents across systematics
 * universes (the distribution that the Pearson correlation coefficient is
 * summarizing in a single number). This provides a check that the bin
 * correlations are adequately captured in a linear approximation.
 *
 * Note that the actual distributions are produced in the `TSCovariance`
 * module; this script only plots them.
 *
 * Usage:
 *
 *   root -l -x 'plot_bincor("input.root")'
 *
 * The input is a ROOT file with covariance matrices, which contains the
 * bin-bin distributions.
 *
 * A. Mastbaum <mastbaum@uchicago.edu>, 2017/10
 */

void plot_bincor(TString infile) {
  TFile* f = TFile::Open(infile);
  assert(f);
  f->ls();

  gStyle->SetOptTitle(0);

  TCanvas* c1 = new TCanvas("c1", "", 1200, 1200);

  // Fake axis labels
  TLatex* l1 = new TLatex(0.05, 0.2, "1#mu1p");
  l1->SetTextAngle(90);
  l1->SetTextFont(132);
  l1->Draw();

  TLatex* l2 = new TLatex(0.2, 0.02, "1#mu1p");
  l2->SetTextFont(132);
  l2->Draw();

  TLatex* l3 = new TLatex(0.05, 0.7, "1e1p");
  l3->SetTextAngle(90);
  l3->SetTextFont(132);
  l3->Draw();

  TLatex* l4 = new TLatex(0.7, 0.02, "1e1p");
  l4->SetTextFont(132);
  l4->Draw();

  TPad* c = new TPad("c", "", 0.07, 0.07, 1.0, 1.0);
  c->Draw();

  // Use a subset of bins (trimming to exclude bins with low/no stats)
  std::vector<size_t> n = {
    2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,0,
    25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46
  };

  c->Divide(n.size(), n.size(), 0, 0);

  char name[50];
  for (long i=0; i<n.size(); i++) {
    for (long j=0; j<n.size(); j++) {
      snprintf(name, 50, "gc_%02lu_%02lu", n[i], n[j]);
      TGraph* gc = (TGraph*) f->Get(name);
      assert(gc);
      TVirtualPad* p = c->cd(n.size()*(n.size()-i-1) + j + 1);
      p->SetBottomMargin(0);
      p->SetTopMargin(0);
      p->SetLeftMargin(0);
      p->SetRightMargin(0);
      gc->SetMarkerStyle(6);
      gc->Draw("ap");
      cout << name << endl;
    }
  }

  c->Update();
}

