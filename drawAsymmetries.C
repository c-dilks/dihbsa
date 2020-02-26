// simply draws all canvases in 'multiGrCanvArr' 
void drawAsymmetries(TString fN = "spinroot/asym_4.root") {
  TFile * f = new TFile(fN,"READ");
  TObjArray * a = (TObjArray*) f->Get("multiGrCanvArr");
  TIter n(a);
  TCanvas * c;
  Int_t i=0;
  TString p = fN;
  p(TRegexp(".root$")) = "";
  TString pn;
  while((c=(TCanvas*)n())) {
    pn = Form("%s.%d.png",p.Data(),i);
    c->Print(pn,"png");
    i++;
  };
};
