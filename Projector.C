// takes spinroot/asym*.root file and rescales the error bars
// according to specifications
//
// arguments:
// - numDays: number of days to project
// - polT: target polarization
// - polE: beam polarization (set to 1 if doing TSA)
// - dilu: dilution factor

Double_t errScale;

void Projector(TString infileN, TString model, Int_t numDays=30,
  Double_t polE = 0.85,
  Double_t polT = 0.85,
  Double_t dilu = 0.2) {

  // input parameters -----------------------
  Double_t evRate = 5.1; // number of events per second
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
  // n.b. pass1 inbending has ~5.3 million dihadrons

  // the statistical uncertainty will be multiplied by this number
  errScale = 1.0 / ( polE * polT * dilu * TMath::Sqrt(evCountProj/evCount) );
  printf("uncertainty will be scaled by: %.2f\n",errScale);


  // models -------------------------------------------
  // - attempt to mimic spectator prediction for G1perp, using
  //   a gaussian envolope * sine function, centered at rho mass
  Float_t sigma = 0.2;
  Float_t Mrho = 0.77;
  TString formuSpec = Form("1/(60*%f)",sigma);
  formuSpec = Form("%s*TMath::Exp(-(x-%f)^2/(2*%f*%f))",
    formuSpec.Data(),Mrho,sigma,sigma);
  formuSpec = Form("%s*TMath::Sin(%f-x)",formuSpec.Data(),Mrho);
  TF1 * funcSpec = new TF1("funcSpec",formuSpec,0,1.5);
  // - A_LU^sinPhiR x-dependence from CLAS12 preliminary data polynomial fit,
  //   used for x-dependence of A_UL^sinPhiR
  Float_t coeff[3] = { 0.00726707, 0.157865, -0.279454 };
  TString formuXdep = Form("%f+%f*x+%f*x*x",coeff[0],coeff[1],coeff[2]);
  TF1 * funcXdep = new TF1("funcXdep",formuXdep,0,1);
  // - flat zero
  TF1 * funcZero = new TF1("funcZero","0",0,1.5);
  // - set model
  TF1 * func;
  if(model=="spec") func = funcSpec;
  else if(model=="xdep") func = funcXdep;
  else func = funcZero;
  // -----------------------------------------------------


  // perform the projection
  TCanvas * canv = new TCanvas("canv","canv",800,600);
  Double_t x,y,ex,ey;
  TLine * zero;
  TObjArrayIter nextGr(asymArr);
  TString grT;
  while(TGraphErrors * gr = (TGraphErrors*) nextGr()) {
    printf("gr = %s = %s\n",gr->GetName(),gr->GetTitle());
    gr->SetMarkerColor(kBlue-3);
    gr->SetLineColor(kBlue-3);
    for(int i=0; i<gr->GetN(); i++) {

      grT = gr->GetTitle();
      grT(TRegexp("LU")) = "UL"; // fix the title LU -> UL
      gr->SetTitle(grT);

      gr->GetPoint(i,x,y);
      ex = gr->GetErrorX(i);
      ey = gr->GetErrorY(i);

      gr->SetPoint(i,x,func->Eval(x));
      gr->SetPointError(i,ex,ey*errScale);

      gr->Draw("APE");
      gr->GetYaxis()->SetRangeUser(-0.05,0.05);
      zero = new TLine(gr->GetXaxis()->GetXmin(),0,
                       gr->GetXaxis()->GetXmax(),0);
      zero->SetLineStyle(2);
      zero->Draw();
      //func->Draw("SAME");
      canv->Print(TString(TString(gr->GetName())+".proj.pdf"),"pdf");

    };
  };
};
