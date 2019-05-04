#include "Asymmetry.h"

ClassImp(Asymmetry)

using namespace std;


Asymmetry::Asymmetry(Int_t phiModulation, Bool_t singleBinMode=false) {
  printf("Instantiating Asymmetry...\n");

  // set up modulation
  if(phiModulation>=0 && phiModulation<nMod) whichMod=phiModulation;
  else {
    fprintf(stderr,"ERROR: bad phiModulation\n");
    return;
  };


  switch(whichMod) {
    case modSinPhiR:
      ModulationTitle = "sin(#phi_{R})";
      modMax = 1;
      break;
    case modSinPhiHR:
      ModulationTitle = "sin(#phi_{h}-#phi_{R})";
      modMax = 1;
      break;
    default:
      fprintf(stderr,"ERROR: bad phiModulation\n");
  };


  // set up (default) binning
  maxIV[vM] = 3;
  maxIV[vX] = 1;
  maxIV[vZ] = 1;
  for(int v=0; v<nIV; v++) nBins[v]=-1;
  // -- mass
  AddBinBound(vM,0);
  AddBinBound(vM,0.5);
  AddBinBound(vM,1);
  AddBinBound(vM,maxIV[vM]);
  // -- x
  AddBinBound(vX,0);
  AddBinBound(vX,0.2);
  AddBinBound(vX,maxIV[vX]);
  // -- z
  AddBinBound(vZ,0);
  AddBinBound(vZ,0.5);
  AddBinBound(vZ,maxIV[vZ]);

  if(singleBinMode) {
    for(int v=0; v<nIV; v++) {
      nBins[v]=-1;
      AddBinBound(v,0);
      AddBinBound(v,maxIV[v]);
    };
    printf("\n-- SINGLE BIN MODE ENABLED\n");
  };

  PrintBinBounds();


  IVname[vM] = "M";
  IVname[vX] = "X";
  IVname[vZ] = "Z";

  IVtitle[vM] = "M_{h}";
  IVtitle[vX] = "x";
  IVtitle[vZ] = "z";

  SpinName[sP] = "P";
  SpinName[sM] = "M";

  SpinTitle[sP] = "spin +";
  SpinTitle[sM] = "spin -";


  ResetVars();

  nEvents = 0;


  // instantiate 1-d distributions
  for(int v=0; v<nIV; v++) {
    for(int b=0; b<nBins[v]; b++) {

      plotName = Form("wDist_%s%d",
       IVname[v].Data(),b
     );
      plotTitle = 
        Form("%s distribution :: %s#in[%.2f,%.2f)",
        IVtitle[v].Data(),
        IVtitle[v].Data(),bound[v][b],bound[v][b+1]
      );
      wDist1[v][b] = new TH1D(plotName,plotTitle,
        w1Bins,bound[v][0],bound[v][nBins[v]]
      );


      for(int s=0; s<nSpin; s++) {

        plotName = Form("mDist%s_%s%d",
          SpinName[s].Data(),
          IVname[v].Data(),b
        );
        plotTitle = 
          Form("%s distribution :: %s :: %s#in[%.2f,%.2f)",
          ModulationTitle.Data(),SpinTitle[s].Data(),
          IVtitle[v].Data(),bound[v][b],bound[v][b+1]
        );
        mDist1[s][v][b] = new TH1D(plotName,plotTitle,
          nModBins,-modMax,modMax);
      };
    };
  };


  // instantiate 2-d distributions
  for(int v1=0; v1<nIV; v1++) {
    for(int v2=0; v2<nIV; v2++) {
      for(int b1=0; b1<nBins[v1]; b1++) {
        for(int b2=0; b2<nBins[v2]; b2++) {

          plotName = Form("wDist_%s%d_%s%d",
            IVname[v1].Data(),b1,
            IVname[v2].Data(),b2
          );
          plotTitle = 
            Form("%s vs %s distribution :: %s#in[%.2f,%.2f)  %s#in[%.2f,%.2f)",
            IVtitle[v2].Data(),IVtitle[v1].Data();
            IVtitle[v1].Data(),bound[v1][b1],bound[v1][b1+1]
            IVtitle[v2].Data(),bound[v2][b2],bound[v2][b2+1]
          );
          wDist2[v1][v2][b1][b2] = new TH1D(plotName,plotTitle,
            w2Bins,bound[v1][0],bound[v1][nBins[v1]],
            w2Bins,bound[v2][0],bound[v2][nBins[v2]]
          );


          for(int s=0; s<nSpin; s++) {

            plotName = Form("mDist%s_%s%d_%s%d",
              SpinName[s].Data(),
              IVname[v1].Data(),b1,
              IVname[v2].Data(),b2
            );
            plotTitle = 
              Form("%s distribution :: %s :: %s#in[%.2f,%.2f)  %s#in[%.2f,%.2f)",
              ModulationTitle.Data(),SpinTitle[s].Data(),
              IVtitle[v1].Data(),bound[v1][b1],bound[v1][b1+1]
              IVtitle[v2].Data(),bound[v2][b2],bound[v2][b2+1]
            );
            mDist2[s][v1][v2][b1][b2] = new TH2D(plotName,plotTitle,
              nModBins,-modMax,modMax);
          };
        }; // eo b2
      }; // eo b1
    }; // eo v2
  }; // eo v1

  // aqui



};



void Asymmetry::AddBinBound(Int_t ivIdx, Float_t newBound) {
  if(ivIdx<0 || ivIdx>=nIV) {
    fprintf(stderr,"ERROR: bad Asymmetry::AddBinBound call");
    return;
  };

  bound[ivIdx][++nBins[ivIdx]] = newBound;

  return;
};


void Asymmetry::PrintBinBounds() {
  for(int v=0; v<nIV; v++) {
    printf("[---] %s bins:  (nbins=%d)\n",IVtitle[v].Data(),nBins[v]);
    for(int b=0; b<nBins[v]; b++) {
      printf(" bin %d:\t\t%.2f\t%.2f\n",b,bound[v][b],bound[v][b+1]);
    };
  };
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
    wDist1[v][binn[v]]->Fill(iv[v]);
    mDist1[spinn][v][binn[v]]->Fill(modulation);
  };

  // fill 2D plots
  /*
  for(int v1=0; v1<nIV; v1++) {
    for(int v2=0; v2<nIV; v2++) {
      wDist2[v1][v2][binn[v1]][binn[v2]]->Fill(iv[v1],iv[v2]);
      mDist2[spinn][v1][v2][binn[v1]][binn[v2]]->Fill(modulation);
    };
  };

  // fill 3D plot
  wDist3[binn[vM]][binn[vX]][binn[vZ]]->Fill(iv[vM],iv[vX],iv[vZ]);
  mDist3[spinn][binn[vM]][binn[vX]][binn[vZ]]->Fill(modulation);
  */


  nEvents++;
};
  
  
    
Float_t Asymmetry::EvalModulation(Float_t PhiH_, Float_t PhiR_) {
  switch(whichMod) {
    case modSinPhiR:
      return TMath::Sin(PhiR_);
      break;
    case modSinPhiHR:
      return TMath::Sin(PhiH_-PhiR_);
      break;
    default:
      fprintf(stderr,"ERROR: bad phiModulation\n");
      return -10000;
  };
};

Int_t Asymmetry::SpinState(Int_t spin_) {
  switch(spin_) {
    case -1: return sM;
    case 1: return sP;
    default:
      fprintf(stderr,"ERROR: bad SpinState request: %d\n",spin_);
      return -10000;
  };
};


void Asymmetry::ResetVars() {
  Mh = -10000;
  x = -10000;
  z = -10000;
  eSpin = 0;
  pSpin = 0;
  PhiH = -10000;
  PhiR = -10000;
  PhPerp = -10000;
  for(int v=0; v<nIV; v++) iv[v]=-10000;
};



void Asymmetry::Write(TFile * f) {
  f->cd();

  // 1D 
  for(int v=0; v<nIV; v++) {
    for(int b=0; b<nBins[v]; b++) {
      wDist1[v][b]->Write();
    };
  };
  for(int v=0; v<nIV; v++) {
    for(int b=0; b<nBins[v]; b++) {
      for(int s=0; s<nSpin; s++) {
        mDist1[s][v][b]->Write();
      };
    };
  };
};




Asymmetry::~Asymmetry() {
};

