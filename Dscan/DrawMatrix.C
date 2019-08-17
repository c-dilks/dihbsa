void DrawMatrix() {
  const Int_t NV = 5;
  TFile * f[NV][NV];
  TString fn;
  TMultiGraph * g[NV][NV];
  TString gn = "multiGr_M";
  TString gt;
  Int_t pad;
  Int_t dsp,dpp;

  TCanvas * canv = new TCanvas("canv","canv",1000,1000);
  canv->Divide(NV,NV);


  for(int v0=0; v0<NV; v0++) {
    for(int v1=0; v1<NV; v1++) {
      dsp = v0 - 2;
      dpp = v1 - 2;
      pad = v0 + NV*v1 + 1;
      canv->cd(pad);
      fn = Form("spin.%d.%d.root",dsp,dpp);
      f[v0][v1] = new TFile(fn,"READ");
      printf("fn=%s @ 0x%p pad=%d\n",fn.Data(),(void*)f[v0][v1],pad);
      if(!(f[v0][v1]->IsZombie())) {
        g[v0][v1] = (TMultiGraph*) f[v0][v1]->Get(gn);
        //gt = TString(g[v0][v1]->GetTitle() + Form("dsp=%d dpp=%d",dsp,dpp));
        g[v0][v1]->SetTitle(Form("dsp=%d dpp=%d",dsp,dpp));
        g[v0][v1]->Draw("LAPE");
      };
    };
  };
};

