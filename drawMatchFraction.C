void drawMatchFraction(TString fN = "match.root") {
  TFile * f = new TFile(fN,"READ");

  TH1D * MhMF[2];
  TH1D * XMF[2];
  TH1D * ZMF[2];
  TString mStr[2];
  mStr[0] = "All";
  mStr[1] = "Matched";
  for(int m=0; m<2; m++) {
    MhMF[m] = (TH1D*) f->Get(TString("MhMF"+mStr[m]));
    XMF[m] = (TH1D*) f->Get(TString("XMF"+mStr[m]));
    ZMF[m] = (TH1D*) f->Get(TString("ZMF"+mStr[m]));
  };

  MhMF[1]->Divide(MhMF[0]);
  XMF[1]->Divide(XMF[0]);
  ZMF[1]->Divide(ZMF[0]);

  new TCanvas();
  MhMF[1]->Draw("E");
  new TCanvas();
  XMF[1]->Draw("E");
  new TCanvas();
  ZMF[1]->Draw("E");
};
