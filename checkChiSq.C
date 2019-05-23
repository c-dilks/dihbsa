void checkChiSq(TString dir="spinout.free") {
  TString fname = dir+"/all.root";
  TFile * infile = new TFile(fname,"READ");
  TH1D * d = (TH1D*)infile->Get("chisqDist");
  Float_t dmin = d->GetXaxis()->GetXmin();
  Float_t dmax = d->GetXaxis()->GetXmax();
  TF1 * f = new TF1("f","[0]*ROOT::Math::chisquared_pdf(x,5)",dmin,dmax);
  new TCanvas();
  d->Rebin(2);
  d->Draw();
  d->Fit(f,"","",dmin,dmax);
};
