#include <cstdlib>
#include <iostream>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TStyle.h"

// DihBsa
#include "Constants.h"
#include "Tools.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"
#include "EventTree.h"
#include "Config.h"


// binning --------------------------------------------
Float_t Q2min = 1e-2; Float_t Q2max = 1e+3;
Float_t xmin = 1e-4;  Float_t xmax = 1;

const Int_t NBINS_Q2 = 4;
Float_t binBounds_Q2[NBINS_Q2] = {
  Q2min,
  2.5,
  10,
  Q2max 
};

const Int_t NBINS_x = 3;
Float_t binBounds_x[NBINS_x] = {
  xmin,
  0.02,
  xmax
};

//-----------------------------------------------------

enum obsEnum {kEle,kDih,kHadA,kHadB,NOBS};
int o;

const Int_t NBINS = NBINS_Q2 * NBINS_x;
Float_t Q2Bin[NBINS][2];
Float_t xBin[NBINS][2];
Double_t numEvents[NBINS];
Double_t maxEvents;
int b;
TString inFiles;
Int_t whichPair;
Int_t whichHad[2];
int h;


TCanvas * PolarCanv(TH2D * dist, TString format);
TCanvas * xQ2Canv(TH2D * dist);
TCanvas * MatrixifyCanv(TCanvas ** canvArr);
TCanvas * MatrixifyDist(TH2F ** distArr, Bool_t logx, Bool_t logy);


int main(int argc, char** argv) {

  // ARGUMENTS
  inFiles = "outroot/*.root";
  whichPair = EncodePairType(kPip,kPim);
  if(argc>1) inFiles = TString(argv[1]);
  if(argc>2) whichPair = (Int_t)strtof(argv[2],NULL);

  // get hadron pair from whichPair; note that in the print out, the 
  // order of hadron 0 and 1 is set by Constants::dihHadIdx
  TString hadTitle[2];
  TString hadName[2];
  TString dihTitle;
  printf("whichPair = 0x%x\n",whichPair);
  DecodePairType(whichPair,whichHad[qA],whichHad[qB]);
  for(h=0; h<2; h++) {
    hadName[h] = PairHadName(whichHad[qA],whichHad[qB],h);
    hadTitle[h] = PairHadTitle(whichHad[qA],whichHad[qB],h);
    printf("hadron %d:  idx=%d  hadron=%s\n",
        h,dihHadIdx(whichHad[qA],whichHad[qB],h),hadTitle[h].Data());
  };
  dihTitle = PairTitle(whichPair);
  printf("dihadron=%s\n",dihTitle.Data());

  // set names and titles for observables
  TString obsN[NOBS];
  TString obsT[NOBS];
  obsN[kEle] = "ele"; obsT[kEle] = "e'"; 
  obsN[kDih] = "dih"; obsT[kDih] = dihTitle;
  obsN[kHadA] = "hadA"; obsT[kHadA] = hadTitle[qA];
  obsN[kHadB] = "hadB"; obsT[kHadB] = hadTitle[qB];


  // setup
  EventTree * ev = new EventTree(inFiles,whichPair);
  Config * conf = new Config();
  gStyle->SetOptStat(0);
  TFile * outfile = new TFile("eicPlots.root","RECREATE");


  // remap all bin bounds into a 1-dimensional array
  // note: the "last" bin for each variable is an extra; it is used to
  //       integrate over other dimensions (see printout)
  for(int bx=0; bx<NBINS_x; bx++) {
    for(int bQ2=0; bQ2<NBINS_Q2; bQ2++) {
      b = bx + bQ2*NBINS_x;
      //b = bx + bQ2*NBINS_x + bPh*NBINS_Q2*NBINS_x; // for another dim
      for(int bb=0; bb<2; bb++) {
        Q2Bin[b][bb] = (bQ2<NBINS_Q2-1) ? binBounds_Q2[bQ2+bb] :
                                          binBounds_Q2[bb==0?0:NBINS_Q2-1];
        xBin[b][bb] = (bx<NBINS_x-1) ? binBounds_x[bx+bb] :
                                          binBounds_x[bb==0?0:NBINS_x-1];
      };
    };
  };
  for(b=0; b<NBINS; b++) {
    printf("b=%d\tQ2:[%g,%g]\tx[%g,%g]\n",
      b,Q2Bin[b][0],Q2Bin[b][1],
      xBin[b][0],xBin[b][1]
    );
  };



  // define plots ---------------------------------------
  const Int_t NPLOTBINS = 50;
  const Int_t N_P_BINS = NPLOTBINS;
  const Int_t N_THETA_BINS = 3*NPLOTBINS;

  // Q2 vs. x planes
  TH2D * Q2vsXfull = new TH2D("Q2vsX_all_events",
    "Log(Q^{2}) vs. Log(x) for all generated events;Log(x);Log(Q^{2})",
    NPLOTBINS,TMath::Log10(xmin),TMath::Log10(xmax),
    NPLOTBINS,TMath::Log10(Q2min),TMath::Log10(Q2max));
  TH2D * Q2vsXcut = new TH2D("Q2vsX_with_cuts",
    "Log(Q^{2}) vs. Log(x) for selected dihadrons;Log(x);Log(Q^{2})",
    NPLOTBINS,TMath::Log10(xmin),TMath::Log10(xmax),
    NPLOTBINS,TMath::Log10(Q2min),TMath::Log10(Q2max));

  // compare (x,Q2) determined from electron to that in pythia event record
  TH2D * Q2vsQ2pythia = new TH2D("Q2vsQ2pythia",
    "Log(Q^{2}) from e' vs. Log(Q^{2}) from pythia event record",
    NPLOTBINS,TMath::Log10(Q2min),TMath::Log10(Q2max),
    NPLOTBINS,TMath::Log10(Q2min),TMath::Log10(Q2max));
  TH2D * XvsXpythia = new TH2D("XvsXpythia",
    "Log(x) from e' vs. Log(x) from pythia event record",
    NPLOTBINS,TMath::Log10(xmin),TMath::Log10(xmax),
    NPLOTBINS,TMath::Log10(xmin),TMath::Log10(xmax));


  // acceptance plots
  TH2D * accPolarLoP[NOBS][NBINS];
  TH2D * accPolarHiP[NOBS][NBINS];
  TH2D * PtVsPzLoP[NOBS][NBINS];
  TH2D * PtVsPzHiP[NOBS][NBINS];
  TH2D * EtaVsP[NOBS][NBINS];
  TH2D * EtaVsPt[NOBS][NBINS];
  TString plotN,plotT,cutT;
  Float_t pMaxLo,pMaxHi;
  for(o=0; o<NOBS; o++) {
    for(b=0; b<NBINS; b++) {

      // (x,Q2) bin range title
      cutT = "";
      if(Q2Bin[b][0]==Q2min && Q2Bin[b][1]==Q2max) cutT += "Full Q^{2}";
      else cutT += Form("%g<Q^{2}<%g",Q2Bin[b][0],Q2Bin[b][1]);
      cutT += " and ";
      if(xBin[b][0]==xmin && xBin[b][1]==xmax) cutT += "Full x";
      else cutT += Form("%g<x<%g",xBin[b][0],xBin[b][1]);

      // bin boundaries
      if(o==kEle) {
        pMaxLo = conf->EbeamEn;
        pMaxHi = 3*pMaxLo;
      } else if(o==kDih || o==kHadA || o==kHadB) {
        pMaxLo = 10;
        pMaxHi = 100;
      };

      // instantiate plots
      plotT = obsT[o] + ", for " + cutT + ";p_{z} [GeV];p_{T} [GeV]";
      plotN = Form("%s_Polar_LoP_%d",obsN[o].Data(),b);
      accPolarLoP[o][b] = new TH2D(plotN,plotT,
        N_THETA_BINS, -PI, PI, N_P_BINS, 0, pMaxLo);
      plotN = Form("%s_Polar_HiP_%d",obsN[o].Data(),b);
      accPolarHiP[o][b] = new TH2D(plotN,plotT,
        N_THETA_BINS, -PI, PI, N_P_BINS, 0, pMaxHi);
      plotN = Form("%s_PtVsPz_LoP_%d",obsN[o].Data(),b);
      PtVsPzLoP[o][b] = new TH2D(plotN,plotT,
        N_P_BINS, -pMaxLo, pMaxLo, N_P_BINS, 0, pMaxLo);
      plotN = Form("%s_PtVsPz_HiP_%d",obsN[o].Data(),b);
      PtVsPzHiP[o][b] = new TH2D(plotN,plotT,
        N_P_BINS, -pMaxHi, pMaxHi, N_P_BINS, 0, pMaxHi);

      plotT = obsT[o] + " #eta vs. p, for " + cutT + ";p [GeV];#eta";
      plotN = Form("%s_EtaVsP_%d",obsN[o].Data(),b);
      EtaVsP[o][b] = new TH2D(plotN,plotT,
        N_P_BINS, 0.01, pMaxLo, N_P_BINS, conf->bdEta[0], conf->bdEta[1]);
      plotT = obsT[o] + " #eta vs. p_{T}, for " + cutT + ";p_{T} [GeV];#eta";
      plotN = Form("%s_EtaVsPt_%d",obsN[o].Data(),b);
      EtaVsPt[o][b] = new TH2D(plotN,plotT,
        N_P_BINS, 0.01, pMaxLo, N_P_BINS, conf->bdEta[0], conf->bdEta[1]);
    };
  };


  // ----------------------------------------------------



  // prepare for event loop
  Bool_t fillPlots;
  Float_t oP,oPt,oTheta,oEta;
  for(b=0; b<NBINS; b++) numEvents[b]=0;


  ///////////////////////////////////////////////
  // EVENT LOOP
  ///////////////////////////////////////////////
  printf("begin loop through %lld events...\n",ev->ENT);
  for(int i=0; i<ev->ENT; i++) {
    //if(i>10000) break; // limiter
    ev->GetEvent(i);

    // fill (x,Q2) plots
    Q2vsXfull->Fill(TMath::Log10(ev->x),TMath::Log10(ev->Q2));
    Q2vsQ2pythia->Fill(TMath::Log10(ev->Q2_pythia),TMath::Log10(ev->Q2));
    XvsXpythia->Fill(TMath::Log10(ev->x_pythia),TMath::Log10(ev->x));
    if(ev->Valid()) {
      Q2vsXcut->Fill(TMath::Log10(ev->x),TMath::Log10(ev->Q2));
    };

    // loop over (x,Q2) bins
    for(b=0; b<NBINS; b++) {

      // override (x,Q2) cuts in EventTree
      ev->cutQ2 = ev->Q2 > Q2Bin[b][0] && ev->Q2 < Q2Bin[b][1];
      ev->cutX = ev->x > xBin[b][0] && ev->x < xBin[b][1];
      ev->cutDIS = ev->cutQ2 && ev->cutX && ev->cutW && ev->cutY;

      fillPlots = ev->Valid(); // use full cuts
      //fillPlots = ev->cutDIS; // use DIS cuts only

      // fill polar plots
      if(fillPlots) {
        for(o=0; o<NOBS; o++) {
          switch(o) {
            case kEle:
              oP = ev->eleP;
              oPt = ev->elePt;
              oTheta = ev->eleTheta;
              oEta = ev->eleEta;
              break;
            case kDih:
              oP = ev->Ph;
              oPt = ev->PhPt;
              oTheta = ev->PhTheta;
              oEta = ev->PhEta;
              break;
            case kHadA:
              oP = ev->hadP[qA];
              oPt = ev->hadPt[qA];
              oTheta = ev->hadTheta[qA];
              oEta = ev->hadEta[qA];
              break;
            case kHadB:
              oP = ev->hadP[qB];
              oPt = ev->hadPt[qB];
              oTheta = ev->hadTheta[qB];
              oEta = ev->hadEta[qB];
              break;
          };

          accPolarLoP[o][b]->Fill(oTheta,oP);
          accPolarHiP[o][b]->Fill(oTheta,oP);
          PtVsPzLoP[o][b]->Fill(oP*TMath::Cos(oTheta),oPt);
          PtVsPzHiP[o][b]->Fill(oP*TMath::Cos(oTheta),oPt);
          EtaVsP[o][b]->Fill(oP,oEta);
          EtaVsPt[o][b]->Fill(oPt,oEta);

        };

        numEvents[b]++;

      };
    };
  };


  // get maximum number of events of all (x,Q2) bins
  maxEvents = 0;
  for(b=0; b<NBINS; b++) {
    printf("numEvents in bin %d: %.f\n",b,numEvents[b]);
    maxEvents = numEvents[b]>maxEvents ? numEvents[b]:maxEvents;
  };
  printf("maxEvents = %.f\n",maxEvents);


  // (x,Q2) plot canvases
  TCanvas * Q2vsXfullCanv = xQ2Canv(Q2vsXfull);
  TCanvas * Q2vsXcutCanv = xQ2Canv(Q2vsXcut);


  // polar plot canvases
  TCanvas * accPolarLoPCanv[NOBS][NBINS];
  TCanvas * accPolarHiPCanv[NOBS][NBINS];
  TCanvas * PtVsPzLoPCanv[NOBS][NBINS];
  TCanvas * PtVsPzHiPCanv[NOBS][NBINS];
  for(o=0; o<NOBS; o++) {
    for(b=0; b<NBINS; b++) {
      accPolarLoPCanv[o][b] = PolarCanv(accPolarLoP[o][b],"pol colz");
      accPolarHiPCanv[o][b] = PolarCanv(accPolarHiP[o][b],"pol colz");
      PtVsPzLoPCanv[o][b] = PolarCanv(PtVsPzLoP[o][b],"colz");
      PtVsPzHiPCanv[o][b] = PolarCanv(PtVsPzHiP[o][b],"colz");
    };
  };


  // matrixify plots
  TCanvas * accPolarLoPMatrix[NOBS];
  TCanvas * accPolarHiPMatrix[NOBS];
  TCanvas * PtVsPzLoPMatrix[NOBS];
  TCanvas * PtVsPzHiPMatrix[NOBS];
  TCanvas * EtaVsPMatrix[NOBS];
  TCanvas * EtaVsPtMatrix[NOBS];
  for(o=0; o<NOBS; o++) {
    accPolarLoPMatrix[o] = MatrixifyCanv(accPolarLoPCanv[o]);
    accPolarHiPMatrix[o] = MatrixifyCanv(accPolarHiPCanv[o]);
    PtVsPzLoPMatrix[o] = MatrixifyCanv(PtVsPzLoPCanv[o]);
    PtVsPzHiPMatrix[o] = MatrixifyCanv(PtVsPzHiPCanv[o]);
    EtaVsPMatrix[o] = MatrixifyDist(EtaVsP[o],1,0);
    EtaVsPtMatrix[o] = MatrixifyDist(EtaVsPt[o],1,0);
  };


  // write output
  Q2vsXfullCanv->Write();
  Q2vsXcutCanv->Write();

  for(o=0; o<NOBS; o++) accPolarLoPMatrix[o]->Write();
  for(o=0; o<NOBS; o++) accPolarHiPMatrix[o]->Write();
  for(o=0; o<NOBS; o++) PtVsPzLoPMatrix[o]->Write();
  for(o=0; o<NOBS; o++) PtVsPzHiPMatrix[o]->Write();
  for(o=0; o<NOBS; o++) EtaVsPMatrix[o]->Write();
  for(o=0; o<NOBS; o++) EtaVsPtMatrix[o]->Write();

  outfile->mkdir("singlePlots");
  outfile->cd("singlePlots");
  Q2vsXfull->Write();
  Q2vsXcut->Write();
  for(o=0; o<NOBS; o++) {
    for(b=0; b<NBINS; b++) accPolarLoPCanv[o][b]->Write();
    for(b=0; b<NBINS; b++) accPolarHiPCanv[o][b]->Write();
    for(b=0; b<NBINS; b++) PtVsPzLoPCanv[o][b]->Write();
    for(b=0; b<NBINS; b++) PtVsPzHiPCanv[o][b]->Write();
    for(b=0; b<NBINS; b++) EtaVsP[o][b]->Write();
    for(b=0; b<NBINS; b++) EtaVsPt[o][b]->Write();
  };
  outfile->cd("/");

  Q2vsQ2pythia->Write();
  XvsXpythia->Write();

  outfile->Close();

};



// take a TH2D and make a polar plot
TCanvas * PolarCanv(TH2D * dist, TString format) {

  // instantiate canvas
  TString canvN = TString(dist->GetName()) + "_Canv";
  TCanvas * canv = new TCanvas(canvN,canvN,800,700);
  canv->SetLogz();
  Double_t radius = dist->GetYaxis()->GetXmax();

  // instantiate frame (because TPad::DrawFrame does not use unique names)
  //canv->DrawFrame(-radius,0,radius,radius,dist->GetTitle());
  TH1F * frame = new TH1F(TString(canvN+"frame"),dist->GetTitle(),
    1000,-radius,radius);
  frame->SetBit(TH1::kNoStats);
  frame->SetMinimum(0);
  frame->SetMaximum(radius);
  frame->GetYaxis()->SetLimits(0,radius);
  frame->Draw(" ");

  // draw circles of constant momentum
  const Int_t nCircles = 3;
  TEllipse * circle[nCircles+1];
  Double_t circleRadius;
  for(int c=0; c<=nCircles; c++) {
    circleRadius = radius - radius*((float)c)/nCircles;
    if(c<nCircles)
      circle[c] = new TEllipse(0,0,circleRadius,circleRadius,0,180);
    else
      circle[c] = new TEllipse(0,0,1.25,1.25,0,180);
    circle[c]->Draw();
  };

  // draw spokes of constant eta
  const Int_t etaSpokeMax = 5;
  TLine * etaSpoke[2*etaSpokeMax+1];
  Double_t theta;
  for(int c=1; c<=etaSpokeMax; c++) {
    theta = Tools::EtaToTheta((Float_t)c);
    etaSpoke[c] = 
     new TLine(0, 0, radius*TMath::Cos(theta), radius*TMath::Sin(theta));
    etaSpoke[c+etaSpokeMax] = 
     new TLine(0, 0, -radius*TMath::Cos(theta), radius*TMath::Sin(theta));
  };
  etaSpoke[0] = new TLine(0, 0, 0, radius);
  for(int c=0; c<2*etaSpokeMax+1; c++) etaSpoke[c]->Draw();

  // draw polar histogram
  dist->SetMinimum(1);
  dist->SetMaximum(maxEvents);
  dist->Draw(TString(format+" same"));
  return canv;
};



// draw bin lines on (x,Q2) plane
TCanvas * xQ2Canv(TH2D * dist) {

  // instantiate canvas
  TString canvN = TString(dist->GetName()) + "_Canv";
  TCanvas * canv = new TCanvas(canvN,canvN,800,700);
  canv->SetLogz();
  dist->Draw("colz");

  // bin lines
  TLine * Q2line[NBINS_Q2];
  TLine * xline[NBINS_x];
  for(b=0; b<NBINS_Q2; b++) {
    Q2line[b] = new TLine(
      TMath::Log10(xmin),
      TMath::Log10(binBounds_Q2[b]),
      TMath::Log10(xmax),
      TMath::Log10(binBounds_Q2[b])
    );
    Q2line[b]->Draw();
  };
  for(b=0; b<NBINS_x; b++) {
    xline[b] = new TLine(
      TMath::Log10(binBounds_x[b]),
      TMath::Log10(Q2min),
      TMath::Log10(binBounds_x[b]),
      TMath::Log10(Q2max)
    );
    xline[b]->Draw();
  };

  return canv;
};


// return a matrix of (x,Q2) bins, containing plots
TCanvas * MatrixifyCanv(TCanvas ** canvArr) {
  TString canvN = TString(canvArr[0]->GetName()) + "Matrix";
  canvN = canvN.ReplaceAll("_0_","_");
  TCanvas * canv = new TCanvas(canvN,canvN,1200,1000);
  canv->Divide(NBINS_x,NBINS_Q2);
  Int_t pad;
  for(b=0; b<NBINS; b++) {
    pad = (NBINS_Q2-b/NBINS_x-1)*NBINS_x + (b%NBINS_x) + 1;
    //printf("b+1=%d pad=%d\n",b+1,pad);
    canv->cd(pad);
    canvArr[b]->DrawClonePad();
  };
  return canv;
};
TCanvas * MatrixifyDist(TH2F ** distArr, Bool_t logx, Bool_t logy) {
  TCanvas * canvases[NBINS];
  TString canvN;
  for(b=0; b<NBINS; b++) {
    canvN = TString(distArr[b]->GetName()) + "Canv";
    canvases[b] = new TCanvas(canvN,canvN,800,700);
    if(logx) canvases[b]->SetLogx();
    if(logy) canvases[b]->SetLogy();
    distArr->Draw("colz");
  };
  return MatrixifyCanv(canvases);
};
