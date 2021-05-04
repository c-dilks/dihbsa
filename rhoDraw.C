void rhoDraw(
 TString infileN = "eicPlots.pythia_5x41_NOsmear.root",
 Int_t binN = 13
) {
  gStyle->SetPalette(kBird);
  TFile *infile = new TFile(infileN,"READ");
  const int N=4;
  TString histN[N] {
    "hadronA_rhoZVsHadZ",
    "hadronA_rhoPperpVsHadPperp",
    "hadronA_rhoQtVsHadQt",
    "hadronA_rhoQtOverQVsHadQtOverQ"
  };
  TString binStr = Form("_%d",binN);
  TCanvas *c;
  TH2D *d;
  TProfile *p;
  for(int i=0; i<N; i++) {
    c = new TCanvas();
    c->SetGrid(1,1);
    d = (TH2D*) infile->Get("singlePlots/"+histN[i]+binStr);
    d->Rebin2D(8,2);
    d->Draw("colz");
    p = d->ProfileX("_pfx",1,-1,"s"); // stddev errors
    p->SetMarkerStyle(kFullCircle);
    p->SetMarkerSize(2);
    p->SetMarkerColor(kBlack);
    p->SetLineColor(kBlack);
    p->SetLineWidth(4);
    p->Draw("same");
  };
  TH1D *m = (TH1D*) infile->Get("singlePlots/dihadron_Mh"+binStr);
  TH1D *mr = (TH1D*) infile->Get("singlePlots/dihadron_rhoMh"+binStr);
  c = new TCanvas();
  m->Draw();
  mr->Draw("same");
  cout << m->GetEntries() << endl;
  cout << mr->GetEntries() << endl;
};
