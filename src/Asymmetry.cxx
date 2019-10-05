#include "Asymmetry.h"

ClassImp(Asymmetry)

using namespace std;


Asymmetry::Asymmetry(
  Binning * binScheme,
  Int_t phiModulation, Int_t dimension, 
  Int_t var0, Int_t bin0,
  Int_t var1, Int_t bin1,
  Int_t var2, Int_t bin2
) {

  // OPTIONS ////////////
  debug = true;
  ///////////////////////


  success = false; // will be true if instance is fully constructed

  // set binning scheme pointer
  BS = binScheme;

  // set up azimuthal modulation
  whichMod = phiModulation;
  modMaxDefault = 1.1;
  asym2d = false;
  switch(whichMod) {
    case modSinPhiR:
      ModulationTitle = "sin(#phi_{R})";
      ModulationName = "sinPhiR";
      modMax = modMaxDefault;
      aziMax = modMaxDefault;
      break;
    case modSinPhiHR:
      ModulationTitle = "sin(#phi_{h}-#phi_{R})";
      ModulationName = "sinPhiHR";
      modMax = modMaxDefault;
      aziMax = modMaxDefault;
      break;
    case weightSinPhiHR:
      ModulationTitle = "p_{T}/M_{h} sin(#phi_{h}-#phi_{R})";
      ModulationName = "sinPhiHRweight";
      modMax = modMaxDefault;
      aziMax = modMaxDefault;
      break;
    case modSinPhiH:
      ModulationTitle = "sin(#phi_{h})";
      ModulationName = "sinPhiH";
      modMax = modMaxDefault;
      aziMax = modMaxDefault;
      break;
    case mod2dSinPhiR:
      asym2d = true;
      ModulationTitle = "2D-fit sin(#phi_{R})";
      ModulationName = "2dPhiHvsPhiR";
      modMax = PI + 0.2;
      aziMax = PI + 0.2;
      break;
    case mod2dWeightSinPhiHR:
      asym2d = true;
      ModulationTitle = "2D-fit p_{T}/M_{h} sin(#phi_{h}-#phi_{R})";
      ModulationName = "2dWeightPhiHvsPhiHR";
      modMax = PI + 0.2;
      aziMax = PI + 0.2;
      break;
    case modTest:
      ModulationTitle = "test modulation"; // for custum testing; search for '+++test'
      ModulationName = "modTest";
      modMax = modMaxDefault;
      aziMax = modMaxDefault;
      break;
    default:
      fprintf(stderr,"ERROR: bad phiModulation\n");
      return;
  };
  if(dimension==-10000) return; // (use this if you only want to do basic things,
                                // like calculate modulations or access modulation
                                // names)

  if(debug) printf("Instantiating Asymmetry...\n");
  printf("  ModulationTitle = %s\n",ModulationTitle.Data());
  printf("  modMax = %f   aziMax = %f\n",modMax,aziMax);

  // set up number of dimensions
  if(dimension>=1 && dimension<=3) whichDim = dimension;
  else {
    fprintf(stderr,"ERROR: bad number of dimensions\n");
    return;
  };



  // check IV mode
  successIVmode = true;
  if(whichDim>=1 && (var0<0 || var0>Binning::nIV) ) successIVmode = false;
  if(whichDim>=2 && (var1<0 || var1>Binning::nIV) ) successIVmode = false;
  if(whichDim>=3 && (var2<0 || var2>Binning::nIV) ) successIVmode = false;
  if(successIVmode) {
    I[0] = var0;  B[0] = bin0;
    I[1] = var1;  B[1] = bin1;
    I[2] = var2;  B[2] = bin2;
  } else {
    fprintf(stderr,"ERROR: bad IV vars\n");
    return;
  };
  if(debug) {
    PrintSettings();
  };



  // set relevant variables for the given IV mode
  for(Int_t d=0; d<3; d++) {
    if(d<whichDim) {
      ivN[d] = BS->IVname[I[d]];
      ivT[d] = BS->IVtitle[I[d]];
      ivMin[d] = BS->minIV[I[d]];
      ivMax[d] = BS->maxIV[I[d]];
    } else {
      ivN[d] = "(unknown)";
      ivT[d] = "(unknown)";
      ivMin[d] = 0;
      ivMax[d] = 0;
    };
  };



  // fixed polarization (for now...)
  pol = 0.86;


  // set up binning title and name
  binT = "::";
  binN = "";
  for(int d=0; d<whichDim; d++) {
    binT += "  " + BS->GetBoundStr(I[d],B[d]);
    binN = Form("%s_%s%d",binN.Data(),ivN[d].Data(),B[d]);
  };
  if(debug) printf("binT = %s\nbinN = %s\n",binT.Data(),binN.Data());
  aName = "A" + binN;


  // ivDist
  // - the histogram of appropriate dimension will be instantiated
  // - the unused dimensions will also be instantiated, but never filled; this is done
  //   to prevent null pointers when streaming to an output file (names prefixed with
  //   "nop"=not operational)
  ivName = Form("ivDist%s",binN.Data());
  if(whichDim == 1) {
    ivTitle = Form("%s distribution %s;%s",
      ivT[0].Data(), binT.Data(),
      ivT[0].Data()
    );
    ivDist1 = new TH1D(ivName,ivTitle,
      iv1Bins,ivMin[0],ivMax[0]
    );
    // nop
    ivDist2 = new TH2D(TString("nop2_"+ivName),TString("nop2"+ivName),1,0,1,1,0,1);
    ivDist3 = new TH3D(TString("nop3_"+ivName),TString("nop3"+ivName),1,0,1,1,0,1,1,0,1);
  }
  else if(whichDim == 2) {
    ivTitle = Form("%s vs. %s %s;%s;%s",
      ivT[1].Data(), ivT[0].Data(), binT.Data(),
      ivT[0].Data(), ivT[1].Data()
    );
    ivDist2 = new TH2D(ivName,ivTitle,
      iv2Bins,ivMin[0],ivMax[0],
      iv2Bins,ivMin[1],ivMax[1]
    );
    // nop
    ivDist1 = new TH1D(TString("nop1_"+ivName),TString("nop1"+ivName),1,0,1);
    ivDist3 = new TH3D(TString("nop3_"+ivName),TString("nop3"+ivName),1,0,1,1,0,1,1,0,1);
  }
  else if(whichDim == 3) {
    ivTitle = Form("%s vs. %s vs. %s %s;%s;%s;%s",
      ivT[2].Data(), ivT[1].Data(), ivT[0].Data(), binT.Data(),
      ivT[0].Data(), ivT[1].Data(), ivT[2].Data()
    );
    ivDist3 = new TH3D(ivName,ivTitle,
      iv3Bins,ivMin[0],ivMax[0],
      iv3Bins,ivMin[1],ivMax[1],
      iv3Bins,ivMin[2],ivMax[2]
    );
    // nop
    ivDist1 = new TH1D(TString("nop1_"+ivName),TString("nop1"+ivName),1,0,1);
    ivDist2 = new TH2D(TString("nop2_"+ivName),TString("nop2"+ivName),1,0,1,1,0,1);
  };
  //if(debug) printf("built %s\n\t%s\n",ivName.Data(),ivTitle.Data());


  // modDist and modBinDist
  modName = Form("modDist%s",binN.Data());
  modTitle = Form("%s distribution %s",ModulationTitle.Data(),binT.Data());
  if(!asym2d) {
    modTitle += ";" + ModulationTitle;
    modDist = new TH1D(modName,modTitle,iv1Bins,-modMax,modMax);
    // nop
    modDist2 = new TH2D(TString("nop_"+modName),TString("nop_"+modName),1,0,1,1,0,1);
  } else {
    modTitle += ";#phi_{R};#phi_{h}";
    modDist2 = new TH2D(modName,modTitle,iv2Bins,-modMax,modMax,iv2Bins,-modMax,modMax);
    // nop
    modDist = new TH1D(TString("nop_"+modName),TString("nop_"+modName),1,0,1);
  };
  //if(debug) printf("built %s\n\t%s\n",modName.Data(),modTitle.Data());

  for(int m=0; m<nModBins; m++) {
    modBinName = Form("%s_bin_%d",modName.Data(),m);
    modBinTitle = Form("bin %d %s",m,modTitle.Data());
    if(!asym2d) {
      modBinDist[m] = new TH1D(modBinName,modBinTitle,iv1Bins,-modMax,modMax);
    } else {
      // nop
      modBinDist[m] = new TH1D(
        TString("nop_"+modBinName),TString("nop_"+modBinName),1,0,1);
    };
  };
  for(int mmH=0; mmH<nModBins2; mmH++) {
    for(int mmR=0; mmR<nModBins2; mmR++) {
      modBinName = Form("%s_bin_H%d_R%d",modName.Data(),mmH,mmR);
      modBinTitle = Form("bin (H%d,R%d) %s",mmH,mmR,modTitle.Data());
      modBinTitle += ";#phi_{R};#phi_{h}";
      if(asym2d) {
        modBinDist2[mmH][mmR] = new TH2D(
          modBinName,modBinTitle,
          iv2Bins,-modMax,modMax,iv2Bins,-modMax,modMax);
      } else {
        //nop
        modBinDist2[mmH][mmR] = new TH2D(
          TString("nop_"+modBinName),TString("nop_"+modBinName),1,0,1,1,0,1);
      };
    };
  };


  // IVvsModDist (only for 1D binning and 1D modulation)
  IVvsModName = Form("IVvsModDist%s",binN.Data());
  IVvsModTitle = Form("%s vs. %s %s;%s;%s",
    ivT[0].Data(),ModulationTitle.Data(),binT.Data(),
    ModulationTitle.Data(),ivT[0].Data()
  );
  if(whichDim==1 && !asym2d) {
    IVvsModDist = new TH2D(IVvsModName,IVvsModTitle,
      iv1Bins,-modMax,modMax,
      iv1Bins,ivMin[0],ivMax[0]
    );
  } else {
    //nop
    IVvsModDist = new TH2D(TString("nop_"+IVvsModName),TString("nop_"+IVvsModName),
      1,0,1,1,0,1);
  };


  // aziDist
  for(int s=0; s<nSpin; s++) {
    aziName[s] = Form("aziDist_%s%s",SpinName(s).Data(),binN.Data());
    aziTitle[s] = Form("%s binned distribution :: %s %s",
      ModulationTitle.Data(),SpinTitle(s).Data(),binT.Data()
    );
    if(!asym2d) {
      aziTitle[s] += ";" + ModulationTitle;
      aziDist[s] = new TH1D(aziName[s],aziTitle[s],nModBins,-aziMax,aziMax);
      //nop
      aziDist2[s] = new TH2D(TString("nop_"+aziName[s]),TString("nop_"+aziName[s]),
        1,0,1,1,0,1);
    } else {
      aziTitle[s] += ";#phi_{R};#phi_{h}"; // PhiR is horizontal, PhiH is vertical
      aziDist2[s] = new TH2D(aziName[s],aziTitle[s],
        nModBins2,-aziMax,aziMax,nModBins2,-aziMax,aziMax);
      //nop
      aziDist[s] = new TH1D(TString("nop_"+aziName[s]),TString("nop_"+aziName[s]),
        1,0,1);
    };
  };

  // yieldDist
  yieldDist = new TH1D(
    TString("yieldDist"+binN),
    TString("yield distribution :: "+binT),
    2,0,2
  );

  // kfDist
  kfLB = 0;
  kfUB = 2;
  kfDist = new TH1D(
    TString("kfDist"+binN),
    TString("K(y) distribution :: "+binT),
    300,kfLB,kfUB);


  // asymGr
  asymName = Form("asym%s",binN.Data());
  if(!asym2d) {
    asymTitle = Form("%s asymmetry %s;%s;asymmetry",
      ModulationTitle.Data(), binT.Data(),  ModulationTitle.Data()
    );
    asymGr = new TGraphErrors();
    asymGr->SetName(asymName);
    asymGr->SetTitle(asymTitle);
    //nop
    asymGr2 = new TGraph2DErrors(); asymGr2->SetTitle(TString("nop_"+asymTitle));
    asymGr2hist = new TH2D(TString("nop_hist"+asymName),TString("nop_hist"+asymName),
      1,0,1,1,0,1);
  } else {
    asymTitle = Form("%s asymmetry %s;#phi_{R};#phi_{h};asymmetry",
      ModulationTitle.Data(), binT.Data()
    );
    asymGr2 = new TGraph2DErrors();
    asymGr2->SetName(asymName);
    asymGr2->SetTitle(asymTitle);
    asymGr2hist = new TH2D(TString("hist"+asymName),asymTitle,
      nModBins2,-aziMax,aziMax,nModBins2,-aziMax,aziMax);
    //nop
    asymGr = new TGraphErrors(); asymGr->SetTitle(TString("nop_"+asymTitle));
  };


  // fit function
  fitFuncName = "fit_"+asymName;
  if(!asym2d) {
    fitFunc = new TF1(fitFuncName,"[0]+[1]*x",-aziMax,aziMax);
    fitFunc->SetParName(0,"B");
    fitFunc->SetParName(1,"A_{LU}");
    //nop
    fitFunc2 = new TF2(TString("nop_"+fitFuncName),"");
  } else {
    if(whichMod == mod2dSinPhiR) {
      fitFunc2 = new TF2(fitFuncName,"[0]*TMath::Sin(x)",
        -aziMax,aziMax,-aziMax,aziMax);
    }
    else if(whichMod == mod2dWeightSinPhiHR) {
      fitFunc2 = new TF2(fitFuncName,"[0]*TMath::Sin(y-x)",
        -aziMax,aziMax,-aziMax,aziMax);
    };
    fitFunc2->SetParName(0,"A_{LU}");
    //nop
    fitFunc = new TF1(TString("nop_"+fitFuncName),"");
  };

  // initialize kinematic variables
  ResetVars();
  nEvents = 0;


  // initialise RooFit
  success = this->InitRooFit();

  if(debug) printf("Asymmetry instantiated.\n");

};



// fill all the plots; to be called in event loop
// -- returns true if filled successfully, or false
//    if not (which usually happens if it's the wrong bin)
Bool_t Asymmetry::AddEvent() {

  // set iv variable
  for(int d=0; d<whichDim; d++) {
    switch(I[d]) {
      case Binning::vM: iv[d] = Mh; break;
      case Binning::vX: iv[d] = x; break;
      case Binning::vZ: iv[d] = z; break;
      case Binning::vPt: iv[d] = PhPerp; break;
      //case Binning::vTh: iv[d] = TMath::Sin(theta); break;
      default: 
        fprintf(stderr,
          "ERROR: Asymmetry::AddEvent does not understand I[%d]=%d\n",d,I[d]);
        return false;
    };
  };


  // check if we're in the proper bin; if not, simply return
  for(int d=0; d<whichDim; d++) {
    if(B[d] != BS->GetBin(I[d],iv[d])) return false;
  };


  // check to make sure kinematics are defined (if they're not, something else
  // probably set them to -10000)
  for(int d=0; d<whichDim; d++) { 
    if(iv[d]<-8000) return KickEvent(TString(ivN[d]+" out of range"),iv[d]);
  };

  if(PhiH<-PIe || PhiH>PIe) return KickEvent("PhiH out of range",PhiH);
  if(PhiR<-PIe || PhiR>PIe) return KickEvent("PhiR out of range",PhiR);

  if(whichMod==modTest && (PhiTest<-PIe || PhiTest>PIe)) 
    return KickEvent("PhiTest out of range",PhiTest);

  if(PhPerp<-8000) return KickEvent("PhPerp out of range",PhPerp);
  if(theta<-0.1 || theta>PIe) return KickEvent("theta out of range",theta);


  // check spin state, which was set by EventTree
  if(spinn<0 || spinn>=nSpin) return false;

  // get kinematic factor
  kf = EvalKinematicFactor();
  if(kf<kfLB || kf>kfUB) return KickEvent("KF out of range",kf);

  // evaluate modulation 
  modulation = EvalModulation(); // (if asym2d==true, modulation=-10000, i.e., not used)

  // set weight (returns 1, unless weighting for G1perp)
  weight = EvalWeight();
  // weight *= kf; // weight events with kinematic factor

  // set RooFit vars
  rfPhiH->setVal(PhiH);
  rfPhiR->setVal(PhiR);
  rfPhiTest->setVal(PhiTest);
  rfWeight->setVal(weight);
  rfTheta->setVal(theta);
  rfSpinCateg->setLabel(rfSpinName[spinn]);
  rfData->add(*rfVars,rfWeight->getVal());






  // fill plots
  // ----------

  if(!asym2d) {
    aziDist[spinn]->Fill(modulation,weight);
    modbin = aziDist[sP]->FindBin(modulation);
    if(modbin>=1 && modbin<=nModBins) {
      modBinDist[modbin-1]->Fill(modulation,weight);
    } else {
      fprintf(stderr,"ERROR: modulation bin not found for modulation=%f\n",modulation);
      return false;
    };
    modDist->Fill(modulation,weight);
  } 
  else {
    aziDist2[spinn]->Fill(PhiR,PhiH,weight);
    modbinR = aziDist2[sP]->GetXaxis()->FindBin(PhiR);
    modbinH = aziDist2[sP]->GetYaxis()->FindBin(PhiH);
    if(modbinR>=1 && modbinR<=nModBins2 && modbinH>=1 && modbinH<=nModBins2) {
      modBinDist2[modbinH-1][modbinR-1]->Fill(PhiR,PhiH,weight);
    } else {
      fprintf(stderr,"ERROR: 2d modulation bin not found (phiH=%f phiR=%f)\n",PhiH,PhiR);
      return false;
    };
    modDist2->Fill(PhiR,PhiH,weight);
  };


  switch(whichDim) {
    case 1: ivDist1->Fill(iv[0]); break;
    case 2: ivDist2->Fill(iv[0],iv[1]); break;
    case 3: ivDist3->Fill(iv[0],iv[1],iv[2]); break;
  };


  if(whichDim==1 && !asym2d) IVvsModDist->Fill(modulation,iv[0],weight);


  yieldDist->Fill(spinn);
  kfDist->Fill(kf);
  
  // increment event counter
  nEvents++;
  return true;
};
  


// calculate the asymmetries; to be called at end of event loop
void Asymmetry::CalculateAsymmetries() {

  
  if(debug) {
    printf("calculate asymmetries for:\n");
    PrintSettings(); 
  };

  // compute relative luminosity
  spinbin = yieldDist->FindBin(sP);
  rNumer = yieldDist->GetBinContent(spinbin);
  spinbin = yieldDist->FindBin(sM);
  rDenom = yieldDist->GetBinContent(spinbin);

  if(rDenom>0) {
    // -- relative luminosity
    rellum = rNumer / rDenom;
    // -- relative luminosity statistical uncertainty (assume poissonian yields)
    rellumErr = sqrt( rNumer * (rNumer+rDenom) / pow(rDenom,3) );
  } else {
    rellum = 0;
    rellumErr = 0;
    fprintf(stderr,"WARNING: rellum denominator==0, abort asym calculation for this bin\n");
    return;
  };
  printf("rellum = %f / %f =  %.3f  +-  %f\n",rNumer,rDenom,rellum,rellumErr);


  // fit asymmetry with RooFit (needs rellum !)
  this->CalculateRooAsymmetries();


   
  // compute asymmetry for each modulation bin
  pointCnt = 0;
  if(!asym2d) {
    for(int m=1; m<=nModBins; m++) {
      yL = aziDist[sP]->GetBinContent(m);
      yR = aziDist[sM]->GetBinContent(m);
      SetAsymGrPoint(m);
    };
  } else {
    for(int mmH=1; mmH<=nModBins2; mmH++) {
      for(int mmR=1; mmR<=nModBins2; mmR++) {
        yL = aziDist2[sP]->GetBinContent(mmR,mmH);
        yR = aziDist2[sM]->GetBinContent(mmR,mmH);
        SetAsymGrPoint(mmH,mmR);
      };
    };
  };


  // fit asymmetry
  if(!asym2d) {
    fitFunc->FixParameter(0,0); // fix fit offset to 0
    asymGr->Fit(fitFunc,"Q","",-aziMax,aziMax);
  } else {
    asymGr2->Fit(fitFunc2,"Q","");
  };

};



// -------------------------------------------
// initialize unbinned ML fit
// -------------------------------------------
Bool_t Asymmetry::InitRooFit() {

  // initialize variables and parameters and data set containers
  // -------------------------------------------------------------
  // - build category to index spins
  rfSpinCateg = new RooCategory("rfSpinCateg","rfSpinCateg");
  for(int s=0; s<nSpin; s++) {
    rfSpinName[s] = "rfSpin" + SpinName(s);
    rfSpinCateg->defineType(rfSpinName[s]);
  };
  // - event vars
  rfPhiH = new RooRealVar("rfPhiH","#phi_{h}",-PIe,PIe);
  rfPhiR = new RooRealVar("rfPhiR","#phi_{R}",-PIe,PIe);
  rfPhiTest = new RooRealVar("rfPhiTest","#phi_{H,sh}",-PIe,PIe);
  rfTheta = new RooRealVar("rfTheta","#theta",-PIe,PIe);
  rfWeight = new RooRealVar("rfWeight","P_{h}^{T}/M_{h}",0,10);

  rfVars = new RooArgSet(*rfPhiH,*rfPhiR,*rfTheta);
  if(whichMod==modTest) rfVars->add(*rfPhiTest);
  rfVars->add(*rfWeight);
  rfVars->add(*rfSpinCateg);

  // - amplitudes (fit parameters)
  rfParamRange = 0.5;
  for(int aa=0; aa<nAmp; aa++) {
    rfAname[aa] = Form("A%d",aa);
    rfA[aa] = new RooRealVar(rfAname[aa],rfAname[aa],-rfParamRange,rfParamRange);
  };
  for(int dd=0; dd<nDparam; dd++) {
    rfDname[dd] = Form("D%d",dd);
    rfD[dd] = new RooRealVar(rfDname[dd],rfDname[dd],-3.0,3.0);
  };
  nAmpUsed = 0;
  nDparamUsed = 0;

  // - yield factor (proportional to actual yield, only for extended fit)
  rfYield[sP] = new RooRealVar("rfYieldP","YP",1e5);
  rfYield[sM] = new RooRealVar("rfYieldM","YM",1e5);
  rfYieldBoth = new RooRealVar("rfYieldBoth","Y",1e5);

  // - polarization and rellum
  rfPol = new RooRealVar("rfPol","P",0,1);
  rfRellum = new RooRealVar("rfRellum","R",0,3);

  // - data sets for each spin
  rfData = new RooDataSet(
    TString("rfData"+binN),TString("rfData"+binN),
    *rfVars,
    rfWeight->GetName()
  );


  // build asymmetry modulation paramaterization "asymExpansion" 
  // = sum { amplitude_i * modulation_i }
  // --------------------------------------------

  // -- modulations
  rfModulation[modSinPhiR] = "TMath::Sin(rfPhiR)";
  rfModulation[modSinPhiHR] = "TMath::Sin(rfPhiH-rfPhiR)";
  rfModulation[weightSinPhiHR] = "TMath::Sin(rfPhiH-rfPhiR)";
  rfModulation[modSinPhiH] = "TMath::Sin(rfPhiH)";
  rfModulation[mod2dSinPhiR] = "TMath::Sin(rfPhiR)";
  rfModulation[mod2dWeightSinPhiHR] = "TMath::Sin(rfPhiH-rfPhiR)";
  rfModulation[modTest] = "TMath::Sin(rfPhiTest)"; // +++test

  // -- define Legendre polynomials; note that these follow the PW expansion of DiFFs,
  //    so overall coefficients may differ slightly from the actual Legendre polynomials
  legendre[0] = "1";
  legendre[1] = "TMath::Cos(rfTheta)";
  legendre[2] = "0.25*(3*TMath::Power(TMath::Cos(rfTheta),2)-1)";

  // -- partial wave expansion factors for asymmetry numerator
  //pwFactorSP = "TMath::Sin(rfTheta)";
  //pwFactorPP = "TMath::Sin(rfTheta)*"+legendre[1];
  pwFactorSP = legendre[0];
  pwFactorPP = legendre[1];

  // -- partial wave expansion factors for asymmetry denominator (unpolarized DiFF)
  pwUnpolDiff = "1+D0*"+legendre[1]+"+D1*"+legendre[2];

  // -- build formula with modulation and PW amplitudes
  Int_t whichFormu = 1; //  <-- <-- <-- <-- <-- <-- <-- <-- <-- <-- <--
  switch(whichFormu) {
    case 0: // test whichMod modulation
      asymExpansion = 
        "A0*" + rfModulation[whichMod];
      rfA[0]->SetTitle("A_{LU}");
      nAmpUsed = 1;
      break;
    case 1: // test linear combination of e(x) and g1perp modulations
      asymExpansion = 
        "A0*" + rfModulation[modSinPhiR] + "+A1*" + rfModulation[weightSinPhiHR];
      rfA[0]->SetTitle("A_{LU}[sin#phi_{R}]");
      rfA[1]->SetTitle("A_{LU}[sin#phi_{hR}]");
      nAmpUsed = 2;
      break;
    case 2: // test single partial wave
      asymExpansion = 
        "A0*" + pwFactorSP + "*" + rfModulation[whichMod];
      rfA[0]->SetTitle("A_{LU}^{sp}[" +ModulationTitle+"]");
      nAmpUsed = 1;
      break;
    case 3: // test 2 partial waves
      asymExpansion = 
        "(A0*" + pwFactorSP + "+A1*" + pwFactorPP + ")*" + rfModulation[whichMod];
      rfA[0]->SetTitle("A_{LU}^{sp}[" +ModulationTitle+"]");
      rfA[1]->SetTitle("A_{LU}^{pp}[" +ModulationTitle+"]");
      nAmpUsed = 2;
      break;
    case 4: // test 2 partial waves with PW-expanded denominator
      asymExpansion = 
        "(A0*" + pwFactorSP + "+A1*" + pwFactorPP + ")*" + rfModulation[whichMod] +
        "/("+pwUnpolDiff+")";
      rfA[0]->SetTitle("A_{LU}^{sp}[" +ModulationTitle+"]");
      rfA[1]->SetTitle("A_{LU}^{pp}[" +ModulationTitle+"]");
      rfD[0]->SetTitle("D^{sp}/D^{ss}");
      rfD[1]->SetTitle("D^{pp}/D^{ss}");

      rfD[0]->setVal(2.0);
      rfD[1]->setVal(-2.0);

      rfD[0]->setConstant(kTRUE);
      rfD[1]->setConstant(kTRUE);

      nAmpUsed = 2;
      nDparamUsed = 2;
      break;
    case 5: // test 2 partial waves and full linear combination of e(x) & g1perp
      asymExpansion = 
        "(A0*" + pwFactorSP + "+A1*" + pwFactorPP + ")*(" +
        "A2*" + rfModulation[modSinPhiR] + "+A3*" + rfModulation[weightSinPhiHR] + ")";
      rfA[0]->SetTitle("a_{LU}^{sp}");
      rfA[1]->SetTitle("a_{LU}^{pp}");
      rfA[2]->SetTitle("a_{LU}[sin#phi_{R}]");
      rfA[3]->SetTitle("a_{LU}[sin#phi_{hR}]");
      nAmpUsed = 4;
      break;
    case 6: // test more linear combination of modulations
      asymExpansion = 
        "A0*" + rfModulation[modSinPhiR] + "+A1*" + rfModulation[weightSinPhiHR] +
        "+A2*TMath::Sin(rfPhiH)+A3*TMath::Sin(rfPhiH+rfPhiR)";
      rfA[0]->SetTitle("A_{LU}[sin#phi_{R}]");
      rfA[1]->SetTitle("A_{LU}[sin(#phi_{h}-#phi_{R}]");
      rfA[2]->SetTitle("A_{LU}[sin#phi_{h}]");
      rfA[3]->SetTitle("A_{LU}[sin(#phi_{h}+#phi_{R}]");
      nAmpUsed = 4;
      break;
    default:
      fprintf(stderr,"ERROR: bad whichFormu\n");
      return false;
  };

  // append polarization factor to asymExpansion
  asymExpansion = "rfPol*("+asymExpansion+")";
      
  // rellum factors
  rellumFactor[sP] = "rfRellum/(rfRellum+1)";
  rellumFactor[sM] = "1/(rfRellum+1)";

  // append yield factor to prefactors
  //for(int s=0; s<nSpin; s++) rellumFactor[s] = "rfYieldBoth*" + rellumFactor[s];
  //rfParams->add(*rfYieldBoth); // DEPRECATED! ADD IT BELOW
  //rellumFactor[sP] = "rfYieldP*" + rellumFactor[sP];
  //rellumFactor[sM] = "rfYieldM*" + rellumFactor[sM];
  //for(int s=0; s<nSpin; s++) rfParams->add(*rfYield[s]); // DEPRECATED! ADD IT BELOW


  // -- build full PDF ( = rellumFactor * ( 1 +/- pol*asymExpansion ) for each spin
  for(int s=0; s<nSpin; s++) {
    rfPdfFormu[s] = rellumFactor[s] + "*(1" + SpinSign(s) + asymExpansion + ")";

    // build list of variables and parameters; we *only* want variables that
    // are actually being used in the PDF
    rfParams[s] = new RooArgSet();
    if(rfPdfFormu[s].Contains("rfPhiH")) rfParams[s]->add(*rfPhiH);
    if(rfPdfFormu[s].Contains("rfPhiR")) rfParams[s]->add(*rfPhiR);
    if(rfPdfFormu[s].Contains("rfPhiTest")) rfParams[s]->add(*rfPhiTest);
    if(rfPdfFormu[s].Contains("rfTheta")) rfParams[s]->add(*rfTheta);
    for(int aa=0; aa<nAmpUsed; aa++) rfParams[s]->add(*rfA[aa]);
    for(int dd=0; dd<nDparamUsed; dd++) rfParams[s]->add(*rfD[dd]);

    rfParams[s]->add(*rfPol);
    rfParams[s]->add(*rfRellum);


    // build pdf
    rfPdf[s] = new RooGenericPdf(
      TString("rfModel" + SpinName(s)),
      TString("rfModel " + SpinTitle(s)),
      rfPdfFormu[s],
      *rfParams[s]
    );
    /*
    rfModelExt[s] = new RooExtendPdf(
      TString("rfPdf" + SpinName(s)),
      TString("rfPdf " + SpinTitle(s)),
      *rfModel[s],
      *rfYield[s]
    );
    rfPdf[s] = new RooAddPdf(
      TString("rfPdf" + SpinName(s)),
      TString("rfPdf " + SpinTitle(s)),
      RooArgList(*rfModel[s])
    );
    */
  };



  // build simultanous PDF 
  rfSimPdf = new RooSimultaneous("rfSimPdf","rfSimPdf",*rfSpinCateg);
  for(int s=0; s<nSpin; s++) rfSimPdf->addPdf(*rfPdf[s],rfSpinName[s]);

  // -log likelihood
  rfNLL = new RooNLLVar("rfNLL","rfNLL",*rfSimPdf,*rfData);
  for(int aa=0; aa<nAmp; aa++) rfNLLplot[aa] = new RooPlot();


  return true;
};


// calculate the asymmetries with RooFit; to be called at end of event loop
void Asymmetry::CalculateRooAsymmetries() {

  // fix polarization and rellum PDF parameters
  rfPol->setVal(pol);
  rfRellum->setVal(rellum);
  rfPol->setConstant(true);
  rfRellum->setConstant(true);


  // fit simultaneous PDF to combined data
  Tools::PrintSeparator(70,"=");
  printf("BEGIN FIT\n");
  Tools::PrintSeparator(70,"=");

  // get number of available threads; if this method fails, set number of threads to 1
  nThreads = (Int_t) std::thread::hardware_concurrency();
  if(nThreads<1) nThreads=1;
  printf("---- fit with %d parallel threads\n",nThreads);

  rfSimPdf->fitTo(*rfData, RooFit::NumCPU(nThreads), RooFit::Save());
  //rfSimPdf->fitTo(*rfData, RooFit::Extended(kTRUE), RooFit::Save(kTRUE));
  //rfSimPdf->fitTo(*rfData,RooFit::PrintLevel(-1));

  // get -log likelihood
  for(int aa=0; aa<nAmpUsed; aa++) {
    rfNLLplot[aa] = rfA[aa]->frame(
      RooFit::Range(-rfParamRange,rfParamRange),
      RooFit::Title(TString("-log(L) scan vs. A"+TString::Itoa(aa,10)))
    );
    rfNLL->plotOn(
      rfNLLplot[aa],
      RooFit::ShiftToZero()
    );
  };


  Tools::PrintSeparator(70,"=");

  // print fit results
  /*
  Tools::PrintTitleBox("ROOFIT RESULTS");
  this->PrintSettings();
  rfResult->Print("v");
  Tools::PrintSeparator(70,"=");
  */

};



// set new asymGr point and error
// -- called by CalculateAsymmetries() for each modulation bin
// -- need to have yL, yR, and rellum set before calling
// -- modBin_ and modBin2_ are used to address modBinDist for getting mean modulation
//    for this modulation bin
void Asymmetry::SetAsymGrPoint(Int_t modBin_, Int_t modBin2_) {

  asymNumer = yL - (rellum * yR);
  asymDenom = yL + (rellum * yR);

  if(asymDenom>0) {
    // compute asymmetry value
    asymVal = (1.0/pol) * (asymNumer/asymDenom);

    // compute asymmetry statistical error
    // -- full formula
    asymErr = ( 2 * rellum * sqrt( yL*pow(yR,2) + yR*pow(yL,2) ) ) / 
              ( pol * pow(yL+rellum*yR,2) );
    // -- compare to simple formula (assumes asym*pol<<1 and R~1)
    //printf("difference = %f\n",asymErr - 1.0 / ( pol * sqrt(yL+yR) ));

    // compute azimuthal modulation value and error
    if(!asym2d) {
      modVal = modBinDist[modBin_-1]->GetMean(); // use modulation bin's mean
      modErr = 0; // for now (TODO)
    } else {
      // using modBinDist2[mmH-1][mmR-1]; x-axis is PhiR; y-axis is PhiH
      modValR = modBinDist2[modBin_-1][modBin2_-1]->GetMean(1);
      modValH = modBinDist2[modBin_-1][modBin2_-1]->GetMean(2);
      modErrR = 0; // for now (TODO)
      modErrH = 0; // for now (TODO)
    };
    

    if(!asym2d) {
      asymGr->SetPoint(pointCnt,modVal,asymVal);
      asymGr->SetPointError(pointCnt,modErr,asymErr);
    } else {
      asymGr2->SetPoint(pointCnt,modValR,modValH,asymVal);
      asymGr2->SetPointError(pointCnt,modErrR,modErrH,asymErr);
      asymGr2hist->SetBinContent(modBin2_,modBin_,asymVal);
    };
    pointCnt++;
  };
};

  
Float_t Asymmetry::EvalModulation() {

  if(PhiH<-1000 || PhiR<-1000) return -10000;

  switch(whichMod) {
    case modSinPhiR:
      return TMath::Sin(PhiR);
      break;
    case modSinPhiHR:
      return TMath::Sin(PhiH-PhiR);
      break;
    case weightSinPhiHR:
      return TMath::Sin(PhiH-PhiR);
      break;
    case modSinPhiH:
      return TMath::Sin(PhiH);
      break;
    case mod2dSinPhiR:
      return -10000; // (not used if asym2d==true)
      break;
    case mod2dWeightSinPhiHR:
      return -10000; // (not used if asym2d==true)
      break;
    case modTest:
      return TMath::Sin(PhiTest); // +++test
      break;
    default:
      fprintf(stderr,"ERROR: bad phiModulation\n");
      return -10000;
  };

};


// g1perp PhPerp/Mh weighting
Float_t Asymmetry::EvalWeight() {
  Float_t wt;
  if( whichMod == weightSinPhiHR
   || whichMod == mod2dWeightSinPhiHR
  ) { 
    wt = Mh>0 ? PhPerp/Mh : 0; 
  }
  else wt = 1;
  return wt;
};


// if e(x) modulation, return W(y)/A(y)
// if G1perp modulation, return C(y)/A(y)
// see EventTree::GetKinematicFactor() for definitions
Float_t Asymmetry::EvalKinematicFactor() {
  switch(whichMod) {
    case modSinPhiR: return kfW / kfA; break;
    case modSinPhiHR: return kfC / kfA; break;
    case weightSinPhiHR: return kfC / kfA; break;
    default: return 1;
  }
};

 
void Asymmetry::ResetVars() {
  Mh = -10000;
  x = -10000;
  z = -10000;
  PhiH = -10000;
  PhiR = -10000;
  PhiTest = -10000;
  PhPerp = -10000;
  theta = -10000;
  spinn = -10000;
  kfA = -10000;
  kfC = -10000;
  kfW = -10000;
  for(int d=0; d<3; d++) iv[d]=-10000;
};


void Asymmetry::PrintSettings() {
  for(Int_t d=0; d<whichDim; d++) printf("  %s bin %d (I[%d]=%d B[%d]=%d)\n",
    (BS->IVname[I[d]]).Data(),B[d],
    d,I[d],d,B[d]
  );
};


// stream pertinent data structures to a TFile
void Asymmetry::StreamData(TFile * tf) {

  tf->cd();
  tf->mkdir(aName);
  tf->cd(TString("/"+aName));

  appName = this->AppFileName(tf);

  printf("writing plots for: "); this->PrintSettings();

  switch(whichDim) {
    case 1: 
      objName = appName + ivDist1->GetName(); ivDist1->Write(objName); 
      break;
    case 2: 
      objName = appName + ivDist2->GetName(); ivDist2->Write(objName); 
      break;
    case 3: 
      objName = appName + ivDist3->GetName(); ivDist3->Write(objName); 
      break;
  };

  if(!asym2d) {
    objName = appName + modDist->GetName(); modDist->Write(objName);
    for(Int_t m=0; m<nModBins; m++) {
      objName = appName + modBinDist[m]->GetName(); modBinDist[m]->Write(objName);
    };
    if(whichDim==1) {
      objName = appName + IVvsModDist->GetName(); IVvsModDist->Write(objName);
    }; 
    for(int s=0; s<nSpin; s++) {
      objName = appName + aziDist[s]->GetName(); aziDist[s]->Write(objName);
    };
  } else {
    objName = appName + modDist2->GetName(); modDist2->Write(objName);
    for(Int_t mmH=0; mmH<nModBins2; mmH++) {
      for(Int_t mmR=0; mmR<nModBins2; mmR++) {
        objName = appName + modBinDist2[mmH][mmR]->GetName();
        modBinDist2[mmH][mmR]->Write(objName);
      };
    };
    for(int s=0; s<nSpin; s++) {
      objName = appName + aziDist2[s]->GetName(); aziDist2[s]->Write(objName);
    };
  };

  objName = appName + yieldDist->GetName(); yieldDist->Write(objName);
  objName = appName + kfDist->GetName(); kfDist->Write(objName);

  objName = appName + rfData->GetName(); rfData->Write(objName);

  tf->cd("/");
};


// append pertinent data structures a TFile to this current instance
void Asymmetry::AppendData(TFile * tf) {
  TH1D * appDist1;
  TH2D * appDist2;
  TH3D * appDist3;
  RooDataSet * appRooDataSet;

  appName = "/" + aName + "/" + this->AppFileName(tf);

  printf("reading plots for: "); this->PrintSettings();

  switch(whichDim) {
    case 1: 
      objName = appName + ivDist1->GetName();
      appDist1 = (TH1D*) tf->Get(objName);
      ivDist1->Add(appDist1); 
      break;
    case 2: 
      objName = appName + ivDist2->GetName();
      appDist2 = (TH2D*) tf->Get(objName);
      ivDist2->Add(appDist2); 
      break;
    case 3: 
      objName = appName + ivDist3->GetName();
      appDist3 = (TH3D*) tf->Get(objName);
      ivDist3->Add(appDist3); 
      break;
  };

  if(!asym2d) {
    objName = appName + modDist->GetName();
    appDist1 = (TH1D*) tf->Get(objName);
    modDist->Add(appDist1);
    for(Int_t m=0; m<nModBins; m++) {
      objName = appName + modBinDist[m]->GetName();
      appDist1 = (TH1D*) tf->Get(objName);
      modBinDist[m]->Add(appDist1);
    };
    if(whichDim==1) {
      objName = appName + IVvsModDist->GetName();
      appDist2 = (TH2D*) tf->Get(objName);
      IVvsModDist->Add(appDist2);
    };
    for(int s=0; s<nSpin; s++) {
      objName = appName + aziDist[s]->GetName();
      appDist1 = (TH1D*) tf->Get(objName);
      aziDist[s]->Add(appDist1);
    };
  } else {
    objName = appName + modDist2->GetName();
    appDist2 = (TH2D*) tf->Get(objName);
    modDist2->Add(appDist2);
    for(Int_t mmH=0; mmH<nModBins2; mmH++) {
      for(Int_t mmR=0; mmR<nModBins2; mmR++) {
        objName = appName + modBinDist2[mmH][mmR]->GetName();
        appDist2 = (TH2D*) tf->Get(objName);
        modBinDist2[mmH][mmR]->Add(appDist2);
      };
    };
    for(int s=0; s<nSpin; s++) {
      objName = appName + aziDist2[s]->GetName();
      appDist2 = (TH2D*) tf->Get(objName);
      aziDist2[s]->Add(appDist2);
    };
  };

  objName = appName + yieldDist->GetName(); 
  appDist1 = (TH1D*) tf->Get(objName);
  yieldDist->Add(appDist1);

  objName = appName + kfDist->GetName(); 
  appDist1 = (TH1D*) tf->Get(objName);
  kfDist->Add(appDist1);

  objName = appName + rfData->GetName();
  appRooDataSet = (RooDataSet*) tf->Get(objName);
  rfData->append(*appRooDataSet); 
  
  tf->cd("/");
};


TString Asymmetry::AppFileName(TFile * tf) {
  TString retstr = TString(tf->GetName());
  retstr(TRegexp("^.*/spin.")) = "stream_";
  retstr(TRegexp(".root$")) = "_";
  return retstr;
};


Bool_t Asymmetry::KickEvent(TString reason,Float_t badValue) {
  fprintf(stderr,"kick event, since %s (value=%f)\n",reason.Data(),badValue);
  return false;
};



Asymmetry::~Asymmetry() {};
