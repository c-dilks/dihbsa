#include "Asymmetry.h"

ClassImp(Asymmetry)

using namespace std;


Asymmetry::Asymmetry(
  Binning * binScheme,
  Int_t phiModulation, Int_t dimension, 
  Int_t var0=0,  Int_t bin0=0,
  Int_t var1=-1, Int_t bin1=-1,
  Int_t var2=-1, Int_t bin2=-1
) {

  printf("Instantiating Asymmetry...\n");
  success = false;
  debug = true;


  // set up number of dimensions
  if(dimensions>=1 && dimensions<=2) whichDim = dimension;
  else {
    fprintf(stderr,"ERROR: bad number of dimensions\n");
    return;
  };


  // set binning scheme pointer
  BS = binScheme;


  // check IV mode
  successIVmode = true;
  if(whichDim>=1 && (var0<0 || var0>Binning:nIV) ) successIVmode = false;
  if(whichDim>=2 && (var1<0 || var1>Binning:nIV) ) successIVmode = false;
  if(whichDim>=3 && (var2<0 || var2>Binning:nIV) ) successIVmode = false;
  if(successIVmode) {
    I[0] = var0;  B[0] = bin0;
    I[1] = var1;  B[1] = bin1;
    I[2] = var2;  B[2] = bin2;
  } else {
    fprintf(stderr,"ERROR: bad IV vars\n");
    return;
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


  // set up azimuthal modulation
  whichMod = phiModulation;
  switch(whichMod) {
    case modSinPhiR:
      ModulationTitle = "sin(#phi_{R})";
      modMax = 1.1;
      break;
    case modSinPhiHR:
      ModulationTitle = "(P_{h}^{perp}/M_{h})sin(#phi_{h}-#phi_{R})";
      modMax = 5;
      break;
    default:
      fprintf(stderr,"ERROR: bad phiModulation\n");
      return;
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
  binT = binT + "  ::";
  if(debug) printf("binT = %s\nbinN = %s\n",binT.Data(),binN.Data());


  // ivDist
  ivName = Form("ivDist%s",binN.Data());
  if(whichDim == 1) {
    ivTitle = Form("%s distribution %s;%s",
      ivT[0].Data(), binT.Data(),
      ivT[0].Data()
    );
    ivDist1 = new TH1D(ivName,ivTitle,
      w1Bins,ivMin[0],ivMax[0]
    );
  }
  else if(whichDim == 2) {
    ivTitle = Form("%s vs. %s %s;%s;%s",
      ivT[1].Data(), ivT[0].Data(), binT.Data(),
      ivT[0].Data(), ivT[1].Data()
    );
    ivDist2 = new TH2D(ivName,ivTitle,
      w2Bins,ivMin[0],ivMax[0],
      w2Bins,ivMin[1],ivMax[1]
    );
  }
  else if(whichDim == 3) {
    ivTitle = Form("%s vs. %s vs. %s %s;%s;%s;%s",
      ivT[2].Data(), ivT[1].Data(), ivT[0].Data(), binT.Data(),
      ivT[0].Data(), ivT[1].Data(), ivT[2].Data()
    );
    ivDist3 = new TH3D(ivName,ivTitle,
      w3Bins,ivMin[0],ivMax[0],
      w3Bins,ivMin[1],ivMax[1],
      w3Bins,ivMin[2],ivMax[2]
    );
  };



  // modDist
  modName = Form("modDist%s",binN.Data());
  modTitle = Form("%s distribution %s;%s",
    ModulationTitle(),binT.Data(),ModulationTitle()
  );
  modDist = new TH1D(modName,modTitle,w1Bins,-modMax,modMax);


  // IVvsModDist (only for 1D!)
  if(whichDim == 1) {
    IVvsModName = Form("IVvsMdist%s",binN.Data());
    IVvsModTitle = Form("%s vs. %s %s;%s;%s",
      ivT[0].Data(),ModulationTitle.Data(),binT.Data(),
      ModulationTitle.Data(),ivT[0]
    );
    IVvsModDist = new TH2D(IVvsModName,IVvsModTitle,
      w1Bins,-modMax,modMax,
      w1Bins,ivMin[0],ivMax[0]
    );
  };


  // aziDist
  for(int s=0; s<nSpin; s++) {
    aziName[s] = Form("aziDist_%s%s",SpinName(s).Data(),binN.Data());
    aziTitle[s] = Form("%s binned distribution :: %s %s;%s",
      ModulationTitle(), SpinTitle(s).Data(), binT.Data(), ModulationTitle()
    );
    aziDist[s] = new TH1D(aziName[s],aziTitle[s],nModBins,-modMax,modMax);
  };


  // asymDist
  asymName = Form("asym%s",binN.Data());
  asymTitle = Form("%s asymmetry %s;%s;asymmetry",
    ModulationTitle.Data(), binT.Data(),  ModulationTitle.Data()
  );
  asymDist = new TGraphErrors();
  asymDist->SetName(asymName);
  asymDist->SetTitle(asymTitle);


  // initialize kinematic variables
  ResetVars();
  nEvents = 0;


  success = true;

};



void Asymmetry::AddBinBound(Int_t ivIdx, Float_t newBound) {
  if(ivIdx<0 || ivIdx>=nIV) {
    fprintf(stderr,"ERROR: bad Asymmetry::AddBinBound call");
    return;
  };

  if(nBins[ivIdx]+1 > nBinsMax) {
    fprintf(stderr,"ERROR: AddBinBound requests more bins than nBinsMax\n");
    return;
  };

  bound[ivIdx][++nBins[ivIdx]] = newBound;

  return;
};


void Asymmetry::PrintBinBounds() {
  printf("\n");
  for(int v=0; v<nIV; v++) {
    printf("[ %s ] %s bins:  (nbins=%d)\n",IVname[v].Data(),IVtitle[v].Data(),nBins[v]);
    for(int b=0; b<nBins[v]; b++) {
      printf(" bin %d:\t\t%.2f\t%.2f\n",b,bound[v][b],bound[v][b+1]);
    };
  };
  printf("\n");
};


Int_t Asymmetry::GetBin(Int_t v_, Float_t iv_) {
  if(iv_<0 || iv_>=nIV) {
    fprintf(stderr,"ERROR: bad Asymmetry::GetBin call\n");
    return -1;
  };

  for(int b=0; b<nBins[v_]; b++) {

    if( iv_ >= bound[v_][b] &&
        iv_ <  bound[v_][b+1] ) {
      return b;
    };
  };

  fprintf(stderr,"ERROR bin not found for iv[%d]=%.2f\n",v_,iv_);
  return -1;
};





void Asymmetry::FillPlots() {
  iv[vM] = Mh;
  iv[vX] = x;
  iv[vZ] = z;
  iv[vPt] = PhPerp;

  // evaluate modulation 
  modulation = EvalModulation(PhiH,PhiR);

  // get spin state number
  spinn = SpinState(eSpin);
  if(spinn<0) return;

  // get bin numbers
  for(int v=0; v<nIV; v++) {
    binn[v] = GetBin(v,iv[v]);
    if(binn[v]<0) return;
  };


  
  // fill 1D plots
  for(int v=0; v<nIV; v++) {
    bDist1[v]->Fill(iv[v]);
    wDist1[v][binn[v]]->Fill(iv[v]);
    mDist1[spinn][v][binn[v]]->Fill(modulation);
  };

  // fill 2D plots
  for(int v1=0; v1<nIV; v1++) {
    for(int v2=0; v2<nIV; v2++) {
      bDist2[v1][v2]->Fill(iv[v1],iv[v2]);
      wDist2[v1][v2][binn[v1]][binn[v2]]->Fill(iv[v1],iv[v2]);
      mDist2[spinn][v1][v2][binn[v1]][binn[v2]]->Fill(modulation);
    };
  };

  // fill 3D plot
  wDist3[binn[vM]][binn[vX]][binn[vZ]]->Fill(iv[vM],iv[vX],iv[vZ]);
  mDist3[spinn][binn[vM]][binn[vX]][binn[vZ]]->Fill(modulation);



  // fill finely-binned modulation dists
  mbDist->Fill(modulation);

  fbin = mDist1[0][0][0]->FindBin(modulation);
  //printf("fbin=%d\n",fbin);
  if(fbin>=1 && fbin<=nModBins) mwDist[fbin-1]->Fill(modulation);
  else {
    fprintf(stderr,
      "ERROR: Asymmetry::FillPlots bad fbin (%d); modulation=%f\n",
      fbin,modulation
    );
  };

  for(int v=0; v<nIV; v++) IVvsMdist[v]->Fill(modulation,iv[v]);


  nEvents++;
};
  


void Asymmetry::CalculateAsymmetries() {
  
  // evaluate 1D asymmetries
  for(int v=0; v<nIV; v++) {
    for(int b=0; b<nBins[v]; b++) {
      EvalAsymmetry(
        asym1[v][b],
        mDist1[sP][v][b],
        mDist1[sM][v][b]
      );
    };
  };

  // evaluate 2D asymmetries
  for(int v1=0; v1<nIV; v1++) {
    for(int v2=0; v2<nIV; v2++) {
      for(int b1=0; b1<nBins[v1]; b1++) {
        for(int b2=0; b2<nBins[v2]; b2++) {
          EvalAsymmetry(
            asym2[v1][v2][b1][b2],
            mDist2[sP][v1][v2][b1][b2],
            mDist2[sM][v1][v2][b1][b2]
          );
        };
      };
    };
  };

  // evaluate 3D asymmetries
  for(int bM=0; bM<nBins[vM]; bM++) {
    for(int bX=0; bX<nBins[vX]; bX++) {
      for(int bZ=0; bZ<nBins[vZ]; bZ++) {
        EvalAsymmetry(
          asym3[bM][bX][bZ],
          mDist3[sP][bM][bX][bZ],
          mDist3[sM][bM][bX][bZ]
        );
      };
    };
  };

};


void Asymmetry::EvalAsymmetry(
  TGraphErrors * asymGr,
  TH1D * mdistL,
  TH1D * mdistR
) {

  if(asymGr==NULL || mdistL==NULL || mdistR==NULL) {
    fprintf(stderr,"ERROR: null pointer in Asymmetry::EvalAsymmetry\n");
    return;
  };


  // compute relative luminosity
  rellumNumer = 0;
  rellumDenom = 0;
  for(int m=1; m<=nModBins; m++) {
    rellumNumer += mdistL->GetBinContent(m);
    rellumDenom += mdistR->GetBinContent(m);
  };
  if(rellumDenom>0) rellum = rellumNumer / rellumDenom;
  else {
    fprintf(stderr,"WARNING: mdistR has 0 yield, abort asym calculation for this bin\n");
    return;
  };

  printf("rellum = %f / %f = %f\n",rellumNumer,rellumDenom,rellum);

   
  // compute asymmetry
  pointCnt = 0;
  for(int m=1; m<=nModBins; m++) {

    yL = mdistL->GetBinContent(m);
    yR = mdistR->GetBinContent(m);

    asymNumer = yL - (rellum * yR);
    asymDenom = yL + (rellum * yR);

    if(asymDenom>0) {
      // compute asymmetry value
      asymVal = (1.0/pol) * (asymNumer/asymDenom);

      // compute asymmetry statistical error
      asymErr = 1.0 / ( pol * TMath::Sqrt(yL+yR) );

      // compute azimuthal modulation value
      //modVal = mdistL->GetBinCenter(m); // use modulation bin's center
      modVal = mwDist[m-1]->GetMean(); // use modulation bin's mean
      
      // compute azimuthal modulation error
      modErr = 0; // azimuthal modulation error

      asymGr->SetPoint(pointCnt,modVal,asymVal);
      asymGr->SetPointError(pointCnt,modErr,asymErr);
      pointCnt++;
    };
  };

  // fit asymmetry
  asymGr->Fit("pol1","Q","",-modMax,modMax);
};


  
    
Float_t Asymmetry::EvalModulation(Float_t PhiH_, Float_t PhiR_) {
  switch(whichMod) {
    case modSinPhiR:
      return TMath::Sin(PhiR_);
      break;
    case modSinPhiHR:
      return (PhPerp/Mh) * TMath::Sin(PhiH_-PhiR_);
      break;
    default:
      fprintf(stderr,"ERROR: bad phiModulation\n");
      return -10000;
  };
};


Int_t Asymmetry::SpinState(Int_t spin_) {
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
  for(int v=0; v<nIV; v++) iv[v]=-10000;
};






Asymmetry::~Asymmetry() {
};

