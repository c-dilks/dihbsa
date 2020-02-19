// prints asymmetry 2d graph, for crosschecking
void grid() {
  TFile * f = new TFile("asym_0.root","READ");
  TGraph2DErrors * g;
  TString gN;
  Double_t phiH,phiR,amp,err;
  for(int i=0; i<2; i++) {
    printf("\nmass bin %d:\n",i);
    gN = Form("asym_M%d",i);
    g = (TGraph2DErrors*) f->Get(gN);
    for(int p=0; p<g->GetN(); p++) {
      g->GetPoint(p,phiR,phiH,amp);
      err = g->GetErrorZ(p);
      printf("%d %f %f %f %f\n",p,phiR,phiH,amp,err);
    };
  };
};
