// takes spinroot/asym*.root file and rescales the error bars
// according to specifications
//
// N.B. the asym*.root file MUST have been produced with polarization==1 and
//      dilution==1
//
// arguments:
// - numDays: number of days to project
// - polT: target polarization
// - polE: beam polarization (set to 1 if doing TSA)
// - dilu: dilution factor

Double_t errScale;

void Projector(TString infileN, TString model,
  TString targetTitle,
  Double_t numDays=60,
  Double_t polE = 0.85,
  Double_t polT = 0.85,
  Double_t dilu = 0.2) {

  // input parameters -----------------------
  Double_t evRate = 5.1; // number of events per second
  Double_t beamEff = 1; // efficiency of data taking (set to 1 for PAC days)
  Bool_t useCustomYields = false; // use custom yields, rather than 
                                  // extrapolating from `infileN`
          if(targetTitle.Contains("He") && model=="specN") 
            useCustomYields = true; // for He3
  // ----------------------------------------
  numDays *= beamEff;

  // get the number of events in asym.root file, and
  // store asymmetry graphs in an array
  TFile * infile = new TFile(infileN,"READ");
  TListIter nextKey(gDirectory->GetListOfKeys());
  TString keyname;
  Double_t evCount = 0;
  TObjArray * asymArr = new TObjArray();
  TString ivBinRange[100];
  int v;
  while(TKey * key = (TKey*) nextKey()) {
    keyname = TString(key->GetName());
    // get number of events
    if(keyname.Contains(TRegexp("^ivFullDist"))) {
      evCount += ((TH1D*)key->ReadObj())->GetEntries();
    }
    // read asymmetry graph
    else if(keyname.Contains(TRegexp("^kindepMA")) &&
           !keyname.Contains("Canv")) {
      asymArr->AddLast((TGraphErrors*)key->ReadObj());
    }
    // get bin ranges
    else if(keyname.Contains(TRegexp("^asym_"))) {
      keyname(TRegexp("^asym_.")) = "";
      v = keyname.Atoi();
      ivBinRange[v] = ((TGraphErrors*)key->ReadObj())->GetTitle();
      ivBinRange[v](TRegexp("^.*::.*in")) = "";
    };
  };

  // number of events expected in numDays
  Double_t evCountProj = numDays * 24 * 60 * 60 * evRate;
  printf("\n");
  printf("projecting %.1f days with %.1f%% efficiency\n",
    numDays/beamEff,beamEff*100);
  printf("number of events in asym.root file: %.0f\n",evCount);
  printf("number of events expected for %.1f days: %.0f\n",numDays,evCountProj);
  printf("yield scale factor: %.2f\n",evCountProj/evCount);
  // n.b. pass1 inbending has ~5.3 million dihadrons



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
  TF1 * funcSpec2 = new TF1("funcSpec2","-0.125*("+formuSpec+")",0,1.5);
  // - A_LU^sinPhiR x-dependence from CLAS12 preliminary data polynomial fit,
  //   used for x-dependence of A_UL^sinPhiR
  Float_t coeff[3] = { 0.00726707, 0.157865, -0.279454 };
  TString formuXdep = Form("%f+%f*x+%f*x*x",coeff[0],coeff[1],coeff[2]);
  TF1 * funcXdep = new TF1("funcXdep",formuXdep,0,1);
  TF1 * funcXdep2 = new TF1("funcXdep2",TString("2*("+formuXdep+")"),0,1);
  // - flat zero
  TF1 * funcZero = new TF1("funcZero","0",0,1.5);
  // - set model
  TF1 * func;
  if(model=="specP") func = funcSpec;
  else if(model=="specN") func = funcSpec2;
  else if(model=="xdepN") func = funcXdep;
  else if(model=="xdepP") func = funcXdep2;
  else func = funcZero;
  // -----------------------------------------------------


  // perform the projection
  TCanvas * canv = new TCanvas("canv","canv",800,600);
  canv->SetGrid(1,1);
  Double_t x,y,ex,ey;
  Double_t ratio;
  TLine * zero;
  TObjArrayIter nextGr(asymArr);
  TString grT;
  TString rootName;
  Bool_t once = true;
  while(TGraphErrors * gr = (TGraphErrors*) nextGr()) {
    //printf("gr = %s = %s\n",gr->GetName(),gr->GetTitle());
    gr->SetMarkerColor(kBlue-3);
    gr->SetLineColor(kBlue-3);
    gr->SetMarkerSize(1.5);
    gr->SetLineWidth(4);

    grT = gr->GetTitle();
    grT(TRegexp("LU")) = "UL"; // fix the title LU -> UL
    if(grT.Contains("theta")) {
      grT(TRegexp("(sin(#theta))\\*(")) = "sin(#theta)";
      grT(TRegexp("(sin(2\\*#theta))\\*(")) = "sin(2*#theta)";
      grT(TRegexp("(pow(sin(#theta),2))\\*(")) = "sin^{2}(#theta)";
      grT(TRegexp("))")) = ")";
    };
    grT = grT + ", for " + targetTitle + " target";
    if(beamEff<1) {
      grT = Form("%s, %d days with %d%% efficiency",
        grT.Data(),(int)(numDays/beamEff),(int)(beamEff*100));
    } else grT = Form("%s, %d PAC days",grT.Data(),(int)(numDays));
    gr->SetTitle(grT);

    for(int i=0; i<gr->GetN(); i++) {

      // the statistical uncertainty will be multiplied by errScale
      if(useCustomYields) {
        // Dien's yield / my yield, for each mass bin, for G1perp projection
        if(x>0.00 && x<0.46)      ratio = 11047542. / 28204.;
        else if(x>0.46 && x<0.60) ratio = 13436566. / 29080.;
        else if(x>0.60 && x<0.72) ratio = 12801074. / 28242.;
        else if(x>0.72 && x<0.81) ratio = 13317744. / 28683.;
        else if(x>0.81 && x<0.93) ratio = 12227050. / 30133.;
        else if(x>0.93 && x<1.10) ratio = 10840390. / 28922.;
        else if(x>1.10 && x<3.00) ratio = 13411398. / 28510.;
        else {
          fprintf(stderr,"ERROR: ratio unknown\n");
          return;
        };
        errScale = 1.0 / ( polE * polT * dilu * TMath::Sqrt(ratio) );
      } else {
        errScale = 1.0 / ( polE * polT * dilu * TMath::Sqrt(evCountProj/evCount) );
      };
      printf("uncertainty will be scaled by: %.2f\n",errScale);
      printf("\n");

      // set projected points
      gr->GetPoint(i,x,y);
      ex = gr->GetErrorX(i);
      ey = gr->GetErrorY(i);
      gr->SetPoint(i,x,func->Eval(x));
      gr->SetPointError(i,ex,ey*errScale);
    };

    // draw projection
    gr->Draw("APE");
    if(model.Contains("xdep")) gr->GetYaxis()->SetRangeUser(-0.01,0.1);
    else gr->GetYaxis()->SetRangeUser(-0.025,0.025);
    zero = new TLine(gr->GetXaxis()->GetXmin(),0,
                     gr->GetXaxis()->GetXmax(),0);
    //zero->SetLineStyle(2);
    zero->SetLineWidth(1);
    zero->Draw();
    //func->Draw("SAME");
    rootName = gr->GetName();
    canv->Print(TString(rootName+".proj.pdf"),"pdf");

    // generate table
    // - bin ranges and means
    if(once) {
      gSystem->RedirectOutput("binranges.tabletex","w");
      for(int i=0; i<gr->GetN(); i++) {
        gr->GetPoint(i,x,y);
        printf("%s & %.2g & %.2g\n",ivBinRange[i].Data(),x,y);
      };
      gSystem->RedirectOutput(0);
      once = false;
      printf("\n");
    };
    // - asymmetry values
    gSystem->RedirectOutput(rootName+".tabletex","w");
    for(int i=0; i<gr->GetN(); i++) {
      gr->GetPoint(i,x,y);
      ey = gr->GetErrorY(i);
      //printf(" & %.2g $\\pm$ %.2g\n",y,ey);
      printf(" & %.2g\n",ey);
    };
    gSystem->RedirectOutput(0);
    printf("\n");
  };
};
