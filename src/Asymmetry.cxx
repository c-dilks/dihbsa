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
  roofitter = true;
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
      ModulationTitle = "weighted sin(#phi_{h}-#phi_{R})";
      ModulationName = "sinPhiHR_weight";
      modMax = modMaxDefault;
      aziMax = modMaxDefault;
      break;
    case modSinPhiH:
      ModulationTitle = "sin(#phi_{h})";
      ModulationName = "sinPhiH";
      modMax = modMaxDefault;
      aziMax = modMaxDefault;
      break;
    case mod2d:
      asym2d = true;
      ModulationTitle = "2D (#phi_{h},#phi_{R})";
      ModulationName = "PhiHvsPhiR";
      modMax = PI + 0.2;
      aziMax = PI + 0.2;
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


  // ivDist
  ivName = Form("ivDist%s",binN.Data());
  if(whichDim == 1) {
    ivTitle = Form("%s distribution %s;%s",
      ivT[0].Data(), binT.Data(),
      ivT[0].Data()
    );
    ivDist1 = new TH1D(ivName,ivTitle,
      iv1Bins,ivMin[0],ivMax[0]
    );
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
  };
  //if(debug) printf("built %s\n\t%s\n",ivName.Data(),ivTitle.Data());


  // modDist and modBinDist
  modName = Form("modDist%s",binN.Data());
  modTitle = Form("%s distribution %s",ModulationTitle.Data(),binT.Data());
  if(!asym2d) {
    modTitle += ";" + ModulationTitle;
    modDist = new TH1D(modName,modTitle,iv1Bins,-modMax,modMax);
  } else {
    modTitle += ";#phi_{R};#phi_{h}";
    modDist2 = new TH2D(modName,modTitle,iv2Bins,-modMax,modMax,iv2Bins,-modMax,modMax);
  };
  //if(debug) printf("built %s\n\t%s\n",modName.Data(),modTitle.Data());

  if(!asym2d) {
    for(int m=0; m<nModBins; m++) {
      modBinName = Form("%s_bin_%d",modName.Data(),m);
      modBinTitle = Form("bin %d %s",m,modTitle.Data());
      modBinTitle += ";" + modTitle;
      modBinDist[m] = new TH1D(
        modBinName.Data(),modBinTitle.Data(),iv1Bins,-modMax,modMax);
    };
  } else {
    for(int mmH=0; mmH<nModBins2; mmH++) {
      for(int mmR=0; mmR<nModBins2; mmR++) {
        modBinName = Form("%s_bin_H%d_R%d",modName.Data(),mmH,mmR);
        modBinTitle = Form("bin (H%d,R%d) %s",mmH,mmR,modTitle.Data());
        modBinTitle += ";#phi_{R};#phi_{h}";
        modBinDist2[mmH][mmR] = new TH2D(
          modBinName.Data(),modBinTitle.Data(),
          iv2Bins,-modMax,modMax,iv2Bins,-modMax,modMax);
      };
    };
  };


  // IVvsModDist (only for 1D binning and 1D modulation)
  if(whichDim==1 && !asym2d) {
    IVvsModName = Form("IVvsModDist%s",binN.Data());
    IVvsModTitle = Form("%s vs. %s %s;%s;%s",
      ivT[0].Data(),ModulationTitle.Data(),binT.Data(),
      ModulationTitle.Data(),ivT[0].Data()
    );
    IVvsModDist = new TH2D(IVvsModName,IVvsModTitle,
      iv1Bins,-modMax,modMax,
      iv1Bins,ivMin[0],ivMax[0]
    );
    //if(debug) printf("built %s\n\t%s\n",IVvsModName.Data(),IVvsModTitle.Data());
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
    } else {
      aziTitle[s] += ";#phi_{R};#phi_{h}"; // PhiR is horizontal, PhiH is vertical
      aziDist2[s] = new TH2D(aziName[s],aziTitle[s],
        nModBins2,-aziMax,aziMax,nModBins2,-aziMax,aziMax);
    };
  };


  // asymGr
  asymName = Form("asym%s",binN.Data());
  if(!asym2d) {
    asymTitle = Form("%s asymmetry %s;%s;asymmetry",
      ModulationTitle.Data(), binT.Data(),  ModulationTitle.Data()
    );
    asymGr = new TGraphErrors();
    asymGr->SetName(asymName);
    asymGr->SetTitle(asymTitle);
  } else {
    asymTitle = Form("%s asymmetry %s;#phi_{R};#phi_{h};asymmetry",
      ModulationTitle.Data(), binT.Data()
    );
    asymGr2 = new TGraph2DErrors();
    asymGr2->SetName(asymName);
    asymGr2->SetTitle(asymTitle);
    asymGr2hist = new TH2D(TString("hist"+asymName),asymTitle,
      nModBins2,-aziMax,aziMax,nModBins2,-aziMax,aziMax);
  };


  // fit function
    fitFuncName = "fit_"+asymName;
  if(!asym2d) {
    fitFunc = new TF1(fitFuncName,"[0]+[1]*x",-aziMax,aziMax);
    fitFunc->SetParName(0,"B");
    fitFunc->SetParName(1,"A_{LU}");
  } else {
    fitFunc2 = new TF2(fitFuncName,"[0]*TMath::Sin(y-x)",-aziMax,aziMax);
    fitFunc2->SetParName(0,"A_{LU}");
  };


  // initialise RooFit
  if(roofitter) {
    // event vars
    rfPhiH = new RooRealVar("rfPhiH","#phi_{h}",-PIe,PIe);
    rfPhiR = new RooRealVar("rfPhiR","#phi_{R}",-PIe,PIe);
    rfWeight = new RooRealVar("rfWeight","P_{h}^{T}/M_{h}",0,4);

    // fit params
    rfAR = new RooRealVar("rfAR","A_{R}",-1,1);

    for(int s=0; s<nSpin; s++) {
      rfData[s] = new RooDataSet(
        TString("rfData"+SpinName(s)),
        TString("rfData"+SpinName(s)),
        RooArgSet(*rfPhiH,*rfPhiR,*rfWeight)
      );
    };
  };



  // initialize kinematic variables
  ResetVars();
  nEvents = 0;
  for(int s=0; s<nSpin; s++) yield[s]=0;

  success = true;

  if(debug) printf("Asymmetry instantiated.\n");

};



// fill all the plots; to be called in event loop
// -- returns true if filled successfully, or false
//    if not (which usually happens if it's the wrong bin)
Bool_t Asymmetry::FillPlots() {

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
          "ERROR: Asymmetry::FillPlots does not understand I[%d]=%d\n",d,I[d]);
        return false;
    };
  };

  // check if we're in the proper bin; if not, simply return
  for(int d=0; d<whichDim; d++) {
    if(B[d] != BS->GetBin(I[d],iv[d])) return false;
  };


  // get spin state number
  spinn = SpinState(eSpin);
  if(spinn<0 || spinn>=nSpin) return false;


  // set RooFit vars
  if(roofitter) {
    *rfPhiH = PhiH;
    *rfPhiR = PhiR;
    *rfWeight = Mh>0 ? PhPerp/Mh : 0;
    rfData[spinn]->add(RooArgSet(*rfPhiH,*rfPhiR,*rfWeight));
  };



  // evaluate modulation 
  modulation = EvalModulation(); // (if asym2d==true, modulation=-10000, i.e., not used)




  // set weight (it's just 1, unless whichMod is set to do a weighted analysis)
  weight = EvalWeight();


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


  
  // increment event counter
  nEvents++;
  yield[spinn]++;
  return true;
};
  


// calculate the asymmetries; to be called at end of event loop
void Asymmetry::CalculateAsymmetries() {

  
  if(debug) {
    printf("calculate asymmetries for:\n");
    PrintSettings(); 
  };

  // compute relative luminosity
  //
  rNumer = yield[sP];
  rDenom = yield[sM];

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
  if(roofitter) this->CalculateRooAsymmetries();


   
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


// calculate the asymmetries with RooFit; to be called at end of event loop
// -- called by CalculateAsymmetries() if roofitter==true
//
void Asymmetry::CalculateRooAsymmetries() {


  // build PDFs for each spin state

  TString rfModulation;
  switch(whichMod) {
    case modSinPhiR:
      rfModulation = "TMath::Sin(rfPhiR)";
      break;
    case modSinPhiHR:
      rfModulation = "TMath::Sin(rfPhiH-rfPhiR)";
      break;
    case weightSinPhiHR:
      rfModulation = "rfWeight*TMath::Sin(rfPhiH-rfPhiR)";
      break;
    case modSinPhiH:
      rfModulation = "TMath::Sin(rfPhiH)";
      break;
    case mod2d:
      // TODO (set to g1perp mod)
      rfModulation = "rfWeight*TMath::Sin(rfPhiH-rfPhiR)";
      break;
    default:
      fprintf(stderr,"ERROR: bad phiModulation\n");
      return;
  };

  rfPdfFormu[sP] = Form("%f*rfAR*%s+1", pol, rfModulation.Data());
  rfPdfFormu[sM] = Form("-(%f/%f)*rfAR*%s+1", pol, rellum, rfModulation.Data());


  rfParams = new RooArgSet(*rfAR,*rfPhiH,*rfPhiR,*rfWeight);


  for(int s=0; s<nSpin; s++) {
    rfPdf[s] = new RooGenericPdf(
      TString("rfPdf"+SpinName(s)),
      TString("rfPdf "+SpinTitle(s)),
      rfPdfFormu[s],
      *rfParams
    );
  };


  // build category and combine data sets from each spin
  rfCateg = new RooCategory("rfCateg","rfCateg");
  for(int s=0; s<nSpin; s++) {
    rfTypeName[s] = "rfType" + SpinName(s);
    rfCateg->defineType(rfTypeName[s]);
  };
  rfCombData = new RooDataSet(
    "rfCombData","rfCombData",
    RooArgSet(*rfPhiH,*rfPhiR,*rfWeight),
    RooFit::Index(*rfCateg),
    RooFit::Import(rfTypeName[sP],*rfData[sP]),
    RooFit::Import(rfTypeName[sM],*rfData[sM])
  );


  // build simultanous PDF 
  rfSimPdf = new RooSimultaneous("rfSimPdf","rfSimPdf",*rfCateg);
  for(int s=0; s<nSpin; s++) rfSimPdf->addPdf(*rfPdf[s],rfTypeName[s]);


  // fit simultaneous PDF to combined data
  rfSimPdf->fitTo(*rfCombData);
  /*
  rfResult = rfSimPdf->fitTo(*rfCombData,RooFit::Save());
  //rfResult = rfSimPdf->fitTo(*rfCombData,RooFit::PrintLevel(-1),RooFit::Save());
  */
  

  // make plot
  for(int s=0; s<nSpin; s++) {
    rfPhiRplot[s] = rfPhiR->frame(
      RooFit::Bins(11),
      RooFit::Title(TString("rfPhiR "+SpinTitle(s)+" Plot"))
    );
    rfCombData->plotOn(
      rfPhiRplot[s],
      RooFit::Cut(TString("rfCateg==rfCateg::"+rfTypeName[s]))
    );
    rfSimPdf->plotOn(
      rfPhiRplot[s],
      RooFit::Slice(*rfCateg,rfTypeName[s]),
      RooFit::ProjWData(*rfCateg,*rfCombData)
    );
  };


  // print fit results
  /*
  Tools::PrintTitleBox("ROOFIT RESULTS");
  this->PrintSettings();
  rfResult->Print();
  //rfResult->Print("v"); // verbose printout
  Tools::PrintSeparator(30);
  */

};



// set new asymGr point and error
// -- called by CalculateAsymmetries() for each modulation bin
// -- need to have yL, yR, and rellum set before calling
// -- modBin_ and modBin2_ are used to address modDistBin for getting mean modulation
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
    case mod2d:
      return -10000; // (not used if asym2d==true)
      break;
    default:
      fprintf(stderr,"ERROR: bad phiModulation\n");
      return -10000;
  };

};


Float_t Asymmetry::EvalWeight() {
  if( whichMod == weightSinPhiHR ||
      whichMod == mod2d
  ) {
    return Mh>0 ? PhPerp/Mh : 0;
  };
  return 1;
};

 
Int_t Asymmetry::SpinState(Int_t spin_) {
  /*
  // +1 -1 convention
  switch(spin_) {
    case 1: return sP;
    case -1: return sM;
    default:
      fprintf(stderr,"WARNING: bad SpinState request: %d\n",spin_);
      return -10000;
  };
  */
  // DNP 2018 convention
  switch(spin_) {
    case 0: return sP;
    case 1: return sM;
    default:
      fprintf(stderr,"WARNING: bad SpinState request: %d\n",spin_);
      return -10000;
  };
};


void Asymmetry::ResetVars() {
  Mh = -10000;
  x = -10000;
  z = -10000;
  eSpin = -10000;
  pSpin = 0;
  PhiH = -10000;
  PhiR = -10000;
  PhPerp = -10000;
  theta = -10000;
  for(int d=0; d<3; d++) iv[d]=-10000;
};


void Asymmetry::PrintSettings() {
  for(Int_t d=0; d<whichDim; d++) printf("  %s bin %d (I[%d]=%d B[%d]=%d)\n",
    (BS->IVname[I[d]]).Data(),B[d],
    d,I[d],d,B[d]
  );
};


Asymmetry::~Asymmetry() {
};

