#include "Asymmetry.h"

ClassImp(Asymmetry)

using namespace std;


Asymmetry::Asymmetry(
  Binning * binScheme,
  Int_t phiModulation, Int_t dimension, 
  Int_t var0=0, Int_t var1=-1, Int_t var2=-1
) {
  printf("Instantiating Asymmetry...\n");
  success = false;


  // set up number of dimensions
  if(dimensions>=1 && dimensions<=2) whichDim = dimension;
  else {
    fprintf(stderr,"ERROR: bad number of dimensions\n");
    return;
  };

  BS = binScheme;

  // check IV mode
  successIVmode = true;
  if(whichDim>=1 && (var0<0 || var0>Binning:nIV) ) successIVmode = false;
  if(whichDim>=2 && (var1<0 || var1>Binning:nIV) ) successIVmode = false;
  if(whichDim>=3 && (var2<0 || var2>Binning:nIV) ) successIVmode = false;
  if(successIVmode) {
    I[0] = var0;
    I[1] = var1;
    I[2] = var2;
  } else {
    fprintf(stderr,"ERROR: bad IV vars\n");
    return;
  };

  // set relevant variables for the given IV mode
  for(Int_t d=0; d<3; d++) {
    nBin[d] = 0;
    ivN[d] = "(unknown)";
    ivT[d] = "(unknown)";
    ivMin[d] = 0;
    ivMax[d] = 1;
  };
  for(Int_t d=0; d<whichDim; d++) {
    nBin[d] = BS->nBins[I[d]];
    ivN[d] = BS->IVname[I[d]];
    ivT[d] = BS->IVtitle[I[d]];
    ivMin[d] = BS->minIV[I[d]];
    ivMax[d] = BS->maxIV[I[d]];
  };

  // set up modulation
  whichMod=phiModulation;
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


  SpinName[sP] = "P";
  SpinName[sM] = "M";

  SpinTitle[sP] = "spin +";
  SpinTitle[sM] = "spin -";


  ResetVars();

  nEvents = 0;


  // instantiate 1-d distributions
  if(whichDim == 1) {

    // ivFullDist
    plotName = Form("ivFullDist_%s",ivN[0].Data());
    plotTitle = Form("%s distribution;%s",ivT[0].Data(),ivT[0].Data());
    ivFullDist1 = new TH1D(plotName,plotTitle,w1Bins,ivMin[0],ivMax[0]);

    // modFullDist
    plotName = Form("modFullDist_%s",ivN[0].Data());
    plotTitle = Form("%s distribution;%s",ModulationTitle.Data(),ModulationTitle.Data());
    modFullDist = new TH1D(plotName,plotTitle,w1Bins,-modMax,modMax);

    // IVvsModDist (only for 1D!)
    plotName = Form("IVvsMdist_%s",ivN[0].Data());
    plotTitle = ivT[0] + " vs. " + ModulationTitle;
    IVvsModDist = new TH2D(plotName,plotTitle,
      w1Bins,-modMax,modMax,
      w1Bins,ivMin[0],ivMax[0]
    );

    for(int b=0; b<nBin[0]; b++) {

      // ivDist
      plotName = Form("ivDist_%s%d",ivN[0].Data(),b);
      plotTitle = Form("%s distribution :: %s;%s",
        ivT[0].Data(),
        (BS->GetBoundStr(I[0],b)).Data(),
        ivT[0].Data()
      );
      ivDist1 = new TH1D(plotName,plotTitle,w1Bins,ivMin[0],ivMax[0]);
      ivVec1.push_back(ivDist1);

      // modDist
      plotName = Form("modDist_%s%d",ivN[0].Data(),b);
      plotTitle = Form("%s distribution :: %s;%s",
        ModulationTitle(),
        (BS->GetBoundStr(I[0],b)).Data(),
        ModulationTitle()
      );
      modDist = new TH1D(plotName,plotTitle,w1Bins,-modMax,modMax);
      modVec.push_back(modDist);

      // asymDist
      plotName = Form("asym_%s%d",ivN[0].Data(),b);
      plotTitle = Form("%s asymmetry :: %s;%s;asym",
        ModulationTitle.Data(),
        (BS->GetBoundStr(I[0],b)).Data(),
        ModulationTitle.Data()
      );
      asymDist = new TGraphErrors();
      asymDist->SetName(plotName);
      asymDist->SetTitle(plotTitle);
      asymVec.push_back(asymDist);

      for(int s=0; s<nSpin; s++) {
        // aziDist
        plotName = Form("aziDist_%s_%s%d",SpinName[s].Data(),ivN[0].Data(),b);
        plotTitle = Form("%s binned distribution :: %s :: %s;%s",
          ModulationTitle(),
          SpinTitle[s].Data(),
          (BS->GetBoundStr(I[0],b)).Data(),
          ModulationTitle()
        );
        aziDist[s] = new TH1D(plotName,plotTitle,nModBins,-modMax,modMax);
        aziVec[s].push_back(aziDist[s]);
      };

    };
  };


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
  /*
  switch(spin_) {
    case -1: return sM;
    case 1: return sP;
    default:
      fprintf(stderr,"WARNING: bad SpinState request: %d\n",spin_);
      return -10000;
  };
  */
  switch(spin_) {
    case 0: return sP; // dnp2019 // TODO
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


void Asymmetry::DrawBoundLines() {

  // 1D 
  for(int v=0; v<nIV; v++) {
    canvName = Form("bCistCanv_%s",IVname[v].Data());
    bDistCanv1[v] = new TCanvas(canvName,canvName,1000,1000);
    bDist1[v]->SetLineColor(kBlack);
    bDist1[v]->SetLineWidth(3);
    bDist1[v]->Draw();

    bMax = bDist1[v]->GetMaximum();
    for(int b=1; b<nBins[v]; b++) {
      boundLine1[v][b] = new TLine(
        bound[v][b], 0,
        bound[v][b], bMax
      );
      boundLine1[v][b]->SetLineColor(kBlue);
      boundLine1[v][b]->SetLineWidth(2);
      boundLine1[v][b]->Draw();
    };

  };


  // 2D 
  for(int v0=0; v0<nIV; v0++) {
    for(int v1=0; v1<nIV; v1++) {
      canvName = Form("bCistCanv_%s_%s",IVname[v0].Data(),IVname[v1].Data());
      bDistCanv2[v0][v1] = new TCanvas(canvName,canvName,1000,1000);
      bDist2[v0][v1]->Draw("colz");

      for(int b=1; b<nBins[v0]; b++) {
        vertLine[v0][v1][b] = new TLine(
          bound[v0][b], bound[v1][0],
          bound[v0][b], bound[v1][nBins[v1]]
        );
        vertLine[v0][v1][b]->SetLineColor(kBlack);
        vertLine[v0][v1][b]->SetLineWidth(3);
        vertLine[v0][v1][b]->Draw();
      };
      
      for(int b=1; b<nBins[v1]; b++) {
        horizLine[v0][v1][b] = new TLine(
          bound[v0][0], bound[v1][b],
          bound[v0][nBins[v0]], bound[v1][b]
        );
        horizLine[v0][v1][b]->SetLineColor(kBlack);
        horizLine[v0][v1][b]->SetLineWidth(3);
        horizLine[v0][v1][b]->Draw();
      };

    };
  };
};


void Asymmetry::WriteObjects(TFile * f) {

  this->DrawBoundLines();

  f->cd();

  // 1D 
  f->mkdir("1D"); f->mkdir("1D/wdists"); f->cd("/1D/wdists");
  for(int v=0; v<nIV; v++) {
    bDistCanv1[v]->Write();
    bDist1[v]->Write();
    for(int b=0; b<nBins[v]; b++) {
      wDist1[v][b]->Write();
    };
  };
  f->cd("/"); f->mkdir("1D/mdists"); f->cd("/1D/mdists");
  for(int v=0; v<nIV; v++) {
    for(int b=0; b<nBins[v]; b++) {
      for(int s=0; s<nSpin; s++) {
        mDist1[s][v][b]->Write();
      };
    };
  };
  f->cd("/"); f->mkdir("1D/asymmetries"); f->cd("/1D/asymmetries");
  for(int v=0; v<nIV; v++) {
    for(int b=0; b<nBins[v]; b++) {
      asym1[v][b]->Write();
    };
  };
  f->cd("/");

  // 2D
  f->mkdir("2D"); f->mkdir("2D/wdists"); f->cd("/2D/wdists");
  for(int v1=0; v1<nIV-1; v1++) {
    for(int v2=v1+1; v2<nIV; v2++) {
      bDistCanv2[v1][v2]->Write();
      bDist2[v1][v2]->Write();
      for(int b1=0; b1<nBins[v1]; b1++) {
        for(int b2=0; b2<nBins[v2]; b2++) {
          wDist2[v1][v2][b1][b2]->Write();
        };
      };
    };
  };
  f->cd("/"); f->mkdir("2D/mdists"); f->cd("/2D/mdists");
  for(int v1=0; v1<nIV-1; v1++) {
    for(int v2=v1+1; v2<nIV; v2++) {
      for(int b1=0; b1<nBins[v1]; b1++) {
        for(int b2=0; b2<nBins[v2]; b2++) {
          for(int s=0; s<nSpin; s++) {
            mDist2[s][v1][v2][b1][b2]->Write();
          };
        };
      };
    };
  };
  f->cd("/"); f->mkdir("2D/asymmetries"); f->cd("/2D/asymmetries");
  for(int v1=0; v1<nIV-1; v1++) {
    for(int v2=v1+1; v2<nIV; v2++) {
      for(int b1=0; b1<nBins[v1]; b1++) {
        for(int b2=0; b2<nBins[v2]; b2++) {
          asym2[v1][v2][b1][b2]->Write();
        };
      };
    };
  };
  f->cd("/");

  // 3D
  f->mkdir("3D"); f->mkdir("3D/wdists"); f->cd("/3D/wdists");
  for(int bM=0; bM<nBins[vM]; bM++) {
    for(int bX=0; bX<nBins[vX]; bX++) {
      for(int bZ=0; bZ<nBins[vZ]; bZ++) {
        wDist3[bM][bX][bZ]->Write();
      };
    };
  };
  f->cd("/"); f->mkdir("3D/mdists"); f->cd("/3D/mdists");
  for(int bM=0; bM<nBins[vM]; bM++) {
    for(int bX=0; bX<nBins[vX]; bX++) {
      for(int bZ=0; bZ<nBins[vZ]; bZ++) {
        for(int s=0; s<nSpin; s++) {
          mDist3[s][bM][bX][bZ]->Write();
        };
      };
    };
  };
  f->cd("/"); f->mkdir("3D/asymmetries"); f->cd("/3D/asymmetries");
  for(int bM=0; bM<nBins[vM]; bM++) {
    for(int bX=0; bX<nBins[vX]; bX++) {
      for(int bZ=0; bZ<nBins[vZ]; bZ++) {
        asym3[bM][bX][bZ]->Write();
      };
    };
  };
  f->cd("/");
 
  mbDist->Write();
  for(int m=0; m<nModBins; m++) mwDist[m]->Write();
  for(int v=0; v<nIV; v++) IVvsMdist[v]->Write();



  printf("\nnEvents = %d\n\n",nEvents);
      
};




Asymmetry::~Asymmetry() {
  printf("destroy Asymmetry instance\n");

  // destroy 1-d distributions
  for(int v=0; v<nIV; v++) {

    if(bDist1[v]) delete bDist1[v];
    if(bDistCanv1[v]) delete bDistCanv1[v];

    for(int b=0; b<nBins[v]; b++) {

      if(wDist1[v][b]) delete wDist1[v][b];
      if(asym1[v][b]) delete asym1[v][b];

      for(int s=0; s<nSpin; s++) {
        if(mDist1[s][v][b]) delete mDist1[s][v][b];
      };

    };
  };


  // destroy 2-d distributions
  for(int v1=0; v1<nIV; v1++) {
    for(int v2=0; v2<nIV; v2++) {

      if(bDist2[v1][v2]) delete bDist2[v1][v2];
      if(bDistCanv2[v1][v2]) delete bDistCanv2[v1][v2];

      for(int b1=0; b1<nBins[v1]; b1++) {
        for(int b2=0; b2<nBins[v2]; b2++) {

          if(wDist2[v1][v2][b1][b2]) delete wDist2[v1][v2][b1][b2];
          if(asym2[v1][v2][b1][b2]) delete asym2[v1][v2][b1][b2];

          for(int s=0; s<nSpin; s++) {
            if(mDist2[s][v1][v2][b1][b2]) delete mDist2[s][v1][v2][b1][b2];
          };

        };
      };
    };
  };


  // destroy 3-d distributions
  for(int bM=0; bM<nBins[vM]; bM++) {
    for(int bX=0; bX<nBins[vX]; bX++) {
      for(int bZ=0; bZ<nBins[vZ]; bZ++) {

        if(wDist3[bM][bX][bZ]) delete wDist3[bM][bX][bZ];
        if(asym3[bM][bX][bZ]) delete asym3[bM][bX][bZ];

        for(int s=0; s<nSpin; s++) {
          if(mDist3[s][bM][bX][bZ]) delete mDist3[s][bM][bX][bZ];
        };

      };
    };
  };


  // destroy lines
  for(int v=0; v<nIV; v++) {
    for(int b=1; b<nBins[v]; b++) {
      if(boundLine1[v][b]) delete boundLine1[v][b];
    };
  };
  for(int v0=0; v0<nIV; v0++) {
    for(int v1=0; v1<nIV; v1++) {
      for(int b=1; b<nBins[v0]; b++) {
        if(vertLine[v0][v1][b]) delete vertLine[v0][v1][b];
      };
      for(int b=1; b<nBins[v1]; b++) {
        if(horizLine[v0][v1][b]) delete horizLine[v0][v1][b];
      };
    };
  };

  // destroy finely-binned modulation dists
  if(mbDist) delete mbDist;
  for(int m=0; m<nModBins; m++) {
    if(mwDist[m]) delete mwDist[m];
  };
  for(int v=0; v<nIV; v++) {
    if(IVvsMdist[v]) delete IVvsMdist[v];
  };
  printf("done\n");

};

