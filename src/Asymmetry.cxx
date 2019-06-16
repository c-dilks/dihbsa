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

  success = false;
  debug = true;


  // set binning scheme pointer
  BS = binScheme;

  // set up azimuthal modulation
  whichMod = phiModulation;
  modMaxDefault = 1.1;
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
    case scaleSinPhiHR:
      ModulationTitle = "(P_{h}^{perp}/M_{h})sin(#phi_{h}-#phi_{R})";
      ModulationName = "sinPhiHR_scale";
      modMax = 5;
      aziMax = BS->GetAziMax(var0,bin0);
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
    ivN[d] = "(unknown)";
    ivT[d] = "(unknown)";
    ivMin[d] = 0;
    ivMax[d] = 0;
  };
  for(Int_t d=0; d<whichDim; d++) {
    ivN[d] = BS->IVname[I[d]];
    ivT[d] = BS->IVtitle[I[d]];
    ivMin[d] = BS->minIV[I[d]];
    ivMax[d] = BS->maxIV[I[d]];
  };




  // fix polarization (for now...)
  pol = 0.86;


  // set up binning title and name
  binT = "::";
  binN = "";
  for(int d=0; d<whichDim; d++) {
    binT = binT + "  " + BS->GetBoundStr(I[d],B[d]);
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
  modTitle = Form("%s distribution %s;%s",
    ModulationTitle.Data(),binT.Data(),ModulationTitle.Data()
  );
  modDist = new TH1D(modName,modTitle,iv1Bins,-modMax,modMax);
  //if(debug) printf("built %s\n\t%s\n",modName.Data(),modTitle.Data());

  for(int m=0; m<nModBins; m++) {
    modBinName[m] = Form("%s_bin%d",modName.Data(),m);
    modBinTitle[m] = Form("bin %d %s",m,modTitle.Data());
    modBinDist[m] = new TH1D(
      modBinName[m].Data(),modBinTitle[m].Data(),iv1Bins,-modMax,modMax);
    //if(debug) printf("built %s\n\t%s\n",modBinName[m].Data(),modBinTitle[m].Data());
  };


  // IVvsModDist (only for 1D!)
  if(whichDim == 1) {
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
    aziTitle[s] = Form("%s binned distribution :: %s %s;%s",
      ModulationTitle.Data(),SpinTitle(s).Data(),binT.Data(),ModulationTitle.Data()
    );
    aziDist[s] = new TH1D(aziName[s],aziTitle[s],nModBins,-aziMax,aziMax);
  };


  // asymGr
  asymName = Form("asym%s",binN.Data());
  asymTitle = Form("%s asymmetry %s;%s;asymmetry",
    ModulationTitle.Data(), binT.Data(),  ModulationTitle.Data()
  );
  asymGr = new TGraphErrors();
  asymGr->SetName(asymName);
  asymGr->SetTitle(asymTitle);

  // fit function
  fitFuncName = "fit_"+asymName;
  fitFunc = new TF1(fitFuncName,"[0]+[1]*x",-aziMax,aziMax);
  fitFunc->SetParName(0,"B");
  fitFunc->SetParName(1,"A_{LU}");


  // initialize kinematic variables
  ResetVars();
  nEvents = 0;
  for(int s=0; s<nSpin; s++) yield[s]=0;

  success = true;

  if(debug) printf("Asymmetry instantiated.\n");

};



// fill all the plots; returns true if filled successfully, or false
// if not (which usually happens if it's the wrong bin)
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


  // evaluate modulation 
  modulation = EvalModulation();


  // get spin state number
  spinn = SpinState(eSpin);
  if(spinn<0) return false;


  // set weight (it's just 1, unless whichMod is set to do a weighted analysis)
  weight = EvalWeight();


  // fill plots

  modbin = aziDist[sP]->FindBin(modulation);
  if(modbin>=1 && modbin<=nModBins) 
    modBinDist[modbin-1]->Fill(modulation,weight);
  else {
    fprintf(stderr,"ERROR: modulation bin not found for modulation=%f\n",modulation);
    return false;
  };

  modDist->Fill(modulation,weight);

  switch(whichDim) {
    case 1: ivDist1->Fill(iv[0]); break;
    case 2: ivDist2->Fill(iv[0],iv[1]); break;
    case 3: ivDist3->Fill(iv[0],iv[1],iv[2]); break;
  };

  aziDist[spinn]->Fill(modulation,weight);

  if(whichDim==1) IVvsModDist->Fill(modulation,iv[0],weight);


  
  // increment event counter
  nEvents++;
  yield[spinn]++;
  return true;
};
  


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
    // -- relative luminosity statistical uncertainty (assume poison yields)
    rellumErr = sqrt( rNumer * (rNumer+rDenom) / pow(rDenom,3) );
  } else {
    rellum = 0;
    rellumErr = 0;
    fprintf(stderr,"WARNING: mdistR has 0 yield, abort asym calculation for this bin\n");
    return;
  };
  printf("rellum = %f / %f =  %.3f  +-  %f\n",rNumer,rDenom,rellum,rellumErr);

   
  // compute asymmetry for each modulation bin
  pointCnt = 0;
  for(int m=1; m<=nModBins; m++) {

    yL = aziDist[sP]->GetBinContent(m);
    yR = aziDist[sM]->GetBinContent(m);

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

      // compute azimuthal modulation value
      modVal = modBinDist[m-1]->GetMean(); // use modulation bin's mean
      
      // compute azimuthal modulation error
      modErr = 0; // azimuthal modulation error

      asymGr->SetPoint(pointCnt,modVal,asymVal);
      asymGr->SetPointError(pointCnt,modErr,asymErr);
      pointCnt++;
    };
  };

  // fit asymmetry
  fitFunc->FixParameter(0,0);
  asymGr->Fit(fitFunc,"Q","",-aziMax,aziMax);
};


  
    
Float_t Asymmetry::EvalModulation() {

  switch(whichMod) {
    case modSinPhiR:
      return TMath::Sin(PhiR);
      break;
    case modSinPhiHR:
      return TMath::Sin(PhiH-PhiR);
      break;
    case scaleSinPhiHR:
      return (PhPerp/Mh) * TMath::Sin(PhiH-PhiR);
      break;
    case weightSinPhiHR:
      return TMath::Sin(PhiH-PhiR);
      break;
    case modSinPhiH:
      return TMath::Sin(PhiH);
      break;
    default:
      fprintf(stderr,"ERROR: bad phiModulation\n");
      return -10000;
  };

};


Float_t Asymmetry::EvalWeight() {
  if(whichMod == weightSinPhiHR) {
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

