// takes spinroot/asym*.root file and rescales the error bars
// according to specifications

Double_t errScale;

void Projector(TString infileN="spinroot/asym_4.root") {

  // input parameters -----------------------
  Int_t numDays = 30; // number of days to extrapolate yields to
  Double_t evRate = 5.1; // number of events per second
  Double_t polT = 0.85; // target polarization
  Double_t polE = 0.85; // beam polarization
  Double_t dilu = 0.2; // dilution factor
  // ----------------------------------------
  
  // get the number of events in asym.root file, and
  // store asymmetry graphs in an array
  TFile * infile = new TFile(infileN,"READ");
  TListIter nextKey(gDirectory->GetListOfKeys());
  TString keyname;
  Double_t evCount = 0;
  TObjArray * asymArr = new TObjArray();
  while(TKey * key = (TKey*) nextKey()) {
    keyname = TString(key->GetName());
    if(keyname.Contains(TRegexp("^ivFullDist"))) {
      evCount += ((TH1D*)key->ReadObj())->GetEntries();
    }
    else if(keyname.Contains(TRegexp("^kindepMA")) &&
           !keyname.Contains("Canv")) {
      asymArr->AddLast((TGraphErrors*)key->ReadObj());
    };
  };

  // number of events expected in numDays
  Double_t evCountProj = numDays * 24 * 60 * 60 * evRate;
  printf("number of events in asym.root file: %.0f\n",evCount);
  printf("number of events expected for %d days: %.0f\n",numDays,evCountProj);
  printf("yield scale factor: %.2f\n",evCountProj/evCount);

  // the statistical uncertainty will be multiplied by this number
  errScale = 1.0 / ( polT * dilu * TMath::Sqrt(evCountProj/evCount) );
  printf("uncertainty will be scaled by: %.2f\n",errScale);

  // perform the projection
  TCanvas * canv = new TCanvas("canv","canv",800,600);
  TObjArrayIter nextGr(asymArr);
  Double_t x,y,ex,ey;
  TLine * zero;
  while(TGraphErrors * gr = (TGraphErrors*) nextGr()) {
    printf("gr = %s = %s\n",gr->GetName(),gr->GetTitle());
    for(int i=0; i<gr->GetN(); i++) {

      gr->GetPoint(i,x,y);
      ex = gr->GetErrorX(i);
      ey = gr->GetErrorY(i);

      gr->SetPoint(i,x,0);
      gr->SetPointError(i,ex,ey);
      //gr->SetPointError(i,ex,ey*errScale);

      gr->Draw("APE");
      gr->GetYaxis()->SetRangeUser(-0.05,0.05);
      zero = new TLine(gr->GetXaxis()->GetXmin(),0,
                       gr->GetXaxis()->GetXmax(),0);
      zero->Draw();
      canv->Print(TString(TString(gr->GetName())+".proj.png"),"png");

    };
  };
};
