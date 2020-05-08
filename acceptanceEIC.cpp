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


// OPTIONS --------------------------------------------
// binning 
Float_t Q2min = 0.1; Float_t Q2max = 2e+3; // Q2min must match generation cut
Float_t xmin = 1e-4;  Float_t xmax = 1;
const Int_t NBINS_Q2 = 3;
const Int_t NBINS_x = 3;
//-----------------------------------------------------
Float_t binBounds_Q2[NBINS_Q2];
Float_t binBounds_x[NBINS_x];


enum obsEnum {kEle,kDih,kHadA,kHadB,NOBS};
int o;

const Int_t NBINS = NBINS_Q2 * NBINS_x;
Float_t Q2Bin[NBINS][2];
Float_t xBin[NBINS][2];
Double_t numEvents[NBINS];
int b;
TString inFiles;
Int_t whichPair;
Int_t whichHad[2];
int h;


TCanvas * PolarCanv(TH2D * dist, TString format, Double_t max);
TCanvas * xQ2Canv(TH2D * dist);
TCanvas * MatrixifyCanv(TCanvas ** canvArr);
TCanvas * MatrixifyDist1(TH1D ** distArr,Bool_t logx,Bool_t logy);
TCanvas * MatrixifyDist2(TH2D ** distArr,Bool_t logx,Bool_t logy,Bool_t logz);
void * FormatDist(TH1 * dist);


int main(int argc, char** argv) {

  // ARGUMENTS
  inFiles = "outroot/*.root";
  whichPair = EncodePairType(kPip,kPim);
  if(argc>1) inFiles = TString(argv[1]);
  if(argc>2) whichPair = (Int_t)strtof(argv[2],NULL);

  // parse beam energies from file name, and define (x,Q2) bins
  // - if reading multiple files, default values listed here will be used
  Int_t EbeamEn = 10;
  Int_t PbeamEn = 100;
  //
  TString filename, token;
  Ssiz_t tf=0;
  if(!(inFiles.Contains("*"))) {
    filename = inFiles;
    printf("filename = %s\n",filename.Data());
    filename(TRegexp("^.*/")) = "";
    filename(TRegexp("\\.root")) = "";
    filename(TRegexp("\\.hepmc")) = "";
    filename(TRegexp("\\.pythia")) = "";
    filename(TRegexp("\\.lund")) = "";
    while(filename.Tokenize(token,tf,"_")) {
      if(token.Contains("x")) {
        sscanf(token.Data(),"%dx%d",&EbeamEn,&PbeamEn);
        break;
      };
    };
    printf("EbeamEn = %d\n",EbeamEn);
    printf("PbeamEn = %d\n",PbeamEn);
  };

  /// define binning
  Float_t Q2mid=5; Float_t xmid=0.02; // 10x100 default
  if( EbeamEn==5 && PbeamEn==41)   { Q2mid=3;    xmid=0.025; };
  if( EbeamEn==5 && PbeamEn==100)  { Q2mid=3;    xmid=0.01;  };
  if( EbeamEn==10 && PbeamEn==100) { Q2mid=5;    xmid=0.02;  };
  if( EbeamEn==18 && PbeamEn==275) { Q2mid=10;   xmid=0.005; };
  binBounds_Q2[0] = Q2min;
  binBounds_Q2[1] = Q2mid;
  binBounds_Q2[2] = Q2max;
  binBounds_x[0] = xmin;
  binBounds_x[1] = xmid;
  binBounds_x[2] = xmax;


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
  obsN[kEle] = "electron"; obsT[kEle] = "e'"; 
  obsN[kDih] = "dihadron"; obsT[kDih] = dihTitle;
  obsN[kHadA] = "hadronA"; obsT[kHadA] = hadTitle[qA];
  obsN[kHadB] = "hadronB"; obsT[kHadB] = hadTitle[qB];


  // setup
  EventTree * ev = new EventTree(inFiles,whichPair);
  Config * conf = new Config();
  gStyle->SetOptStat(0);
  TString outfileN= inFiles;
  if(inFiles.Contains("*")) outfileN = "eicPlots.root";
  else outfileN(TRegexp("^.*/")) = "eicPlots.";
  printf("outfileN = %s\n",outfileN.Data());
  TFile * outfile = new TFile(outfileN,"RECREATE");


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
    "Q^{2} vs. x for all generated events;x;Q^{2} [GeV^{2}]",
    NPLOTBINS,xmin,xmax,NPLOTBINS,Q2min,Q2max);
  TH2D * Q2vsXcut = new TH2D("Q2vsX_with_cuts",
    "Q^{2} vs. x for selected dihadrons;x;Q^{2} [GeV^{2}]",
    NPLOTBINS,xmin,xmax,NPLOTBINS,Q2min,Q2max);
  Tools::BinLog(Q2vsXfull->GetXaxis());
  Tools::BinLog(Q2vsXfull->GetYaxis());
  Tools::BinLog(Q2vsXcut->GetXaxis());
  Tools::BinLog(Q2vsXcut->GetYaxis());

  // compare (x,Q2) determined from electron to that in pythia event record
  TH2D * Q2vsQ2pythia = new TH2D("Q2vsQ2pythia",
    "Q^{2} from e' vs. Q^{2} from pythia event record",
    NPLOTBINS,Q2min,Q2max,NPLOTBINS,Q2min,Q2max);
  TH2D * XvsXpythia = new TH2D("XvsXpythia",
    "x from e' vs. x from pythia event record",
    NPLOTBINS,xmin,xmax,NPLOTBINS,xmin,xmax);
  Tools::BinLog(Q2vsQ2pythia->GetXaxis());
  Tools::BinLog(Q2vsQ2pythia->GetYaxis());
  Tools::BinLog(XvsXpythia->GetXaxis());
  Tools::BinLog(XvsXpythia->GetYaxis());

  // acceptance plots
  TH2D * accPolarLoP[NOBS][NBINS];
  TH2D * accPolarHiP[NOBS][NBINS];
  TH2D * PtVsPzLoP[NOBS][NBINS];
  TH2D * PtVsPzHiP[NOBS][NBINS];
  TH2D * EtaVsP[NOBS][NBINS];
  TH2D * EtaVsPt[NOBS][NBINS];
  TH2D * EtaVsZ[NOBS][NBINS];
  TH2D * PhPerpVsPt[NOBS][NBINS]; // dihadron PhPerp vs. obs pT
  TH2D * QtVsPt[NOBS][NBINS]; // qT = PhPerp/z, vs. obs pT
  TH2D * QtVsEta[NOBS][NBINS]; // qT = PhPerp/z, vs. obs pT
  TH2D * PtVsY[NOBS][NBINS];
  // - dihadrons only
  TH2D * PhiHvsPhiR[NBINS];
  TH2D * corrXF[NBINS];
  TH1D * distMh[NBINS];
  TH1D * distMx[NBINS];
  TH1D * distTheta[NBINS];
  TH2D * PhiHvsP[NBINS];
  TH2D * PhiRvsP[NBINS];
  TH2D * thetaVsP[NBINS];
  TH2D * QtVsY[NBINS];
  // - electrons only
  TH1D * distY[NBINS];
  TH1D * distW[NBINS];
  ////
  TString plotN,plotT,cutT;
  Float_t pMaxLo,pMaxHi;
  Float_t ptMin,ptMax;
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
        pMaxLo = (Float_t) EbeamEn;
        pMaxHi = 3*pMaxLo;
        ptMin = 0.1;
        ptMax = 1.5*pMaxLo;
      } else if(o==kDih || o==kHadA || o==kHadB) {
        pMaxLo = 10;
        pMaxHi = 50;
        ptMin = 0.001;
        ptMax = pMaxLo;
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
        N_P_BINS, 0.5, 3*pMaxLo,
        N_P_BINS, conf->bdEta[0], conf->bdEta[1]);
      Tools::BinLog(EtaVsP[o][b]->GetXaxis());

      plotT = obsT[o] + " #eta vs. p_{T,lab}, for " + cutT + 
        ";p_{T,lab} [GeV];#eta";
      plotN = Form("%s_EtaVsPt_%d",obsN[o].Data(),b);
      EtaVsPt[o][b] = new TH2D(plotN,plotT,
        N_P_BINS, ptMin, ptMax,
        N_P_BINS, conf->bdEta[0], conf->bdEta[1]);
      Tools::BinLog(EtaVsPt[o][b]->GetXaxis());

      plotT = obsT[o] + " #eta vs. z, for " + cutT + ";z;#eta";
      plotN = Form("%s_EtaVsZ_%d",obsN[o].Data(),b);
      EtaVsZ[o][b] = new TH2D(plotN,plotT,
        N_P_BINS, 0, 1, N_P_BINS, conf->bdEta[0], conf->bdEta[1]);

      plotT = obsT[kDih] + " P_{h}^{perp} vs. " + obsT[o] + 
        " p_{T,lab}, for " + cutT + 
        ";" + obsT[o] + " p_{T,lab} [GeV];P_{h}^{perp} [GeV]";
      plotN = Form("%s_PhPerpVsPt_%d",obsN[o].Data(),b);
      PhPerpVsPt[o][b] = new TH2D(plotN,plotT,
        N_P_BINS, ptMin, ptMax,
        N_P_BINS, ptMin, 10*ptMax);
      Tools::BinLog(PhPerpVsPt[o][b]->GetXaxis());
      Tools::BinLog(PhPerpVsPt[o][b]->GetYaxis());

      plotT = "q_{T} vs. " + obsT[o] +
        " p_{T,lab}, for " + cutT + 
        ";" + obsT[o] + " p_{T,lab} [GeV];q_{T} [GeV]";
      plotN = Form("%s_QtVsPt_%d",obsN[o].Data(),b);
      QtVsPt[o][b] = new TH2D(plotN,plotT,
        N_P_BINS, ptMin, ptMax,
        N_P_BINS, ptMin, 100*ptMax);
      Tools::BinLog(QtVsPt[o][b]->GetXaxis());
      Tools::BinLog(QtVsPt[o][b]->GetYaxis());

      plotT = "q_{T} vs. " + obsT[o] +
        " #eta, for " + cutT + 
        ";" + obsT[o] + " #eta;q_{T} [GeV]";
      plotN = Form("%s_QtVsEta_%d",obsN[o].Data(),b);
      QtVsEta[o][b] = new TH2D(plotN,plotT,
        N_P_BINS, conf->bdEta[0], conf->bdEta[1],
        N_P_BINS, ptMin, 100*ptMax);
      Tools::BinLog(QtVsEta[o][b]->GetYaxis());

      plotT = "p_{T,lab} vs. " + obsT[o] +
        " y, for " + cutT + 
        ";" + obsT[o] + " y;p_{T,lab} [GeV]";
      plotN = Form("%s_PtVsY_%d",obsN[o].Data(),b);
      PtVsY[o][b] = new TH2D(plotN,plotT,
        N_P_BINS, 1e-3, 1,
        N_P_BINS, ptMin, ptMax);
      Tools::BinLog(PtVsY[o][b]->GetXaxis());
      Tools::BinLog(PtVsY[o][b]->GetYaxis());

      if(o==kDih) {
        plotT = obsT[o] + " #phi_{h} vs. #phi_{R}, for " + cutT + 
          ";#phi_{R};#phi_{h}";
        plotN = Form("%s_PhiHvsPhiR_%d",obsN[o].Data(),b);
        PhiHvsPhiR[b] = new TH2D(plotN,plotT,
          N_P_BINS, -PI, PI, N_P_BINS, -PI, PI);

        plotT = hadTitle[qA] + " x_{F} vs. " + 
                hadTitle[qB] + " x_{F}, for " + cutT +
                ";" + hadTitle[qB] + " x_{F}" +
                ";" + hadTitle[qA] + " x_{F}";
        plotN = Form("%s_corrXF_%d",obsN[o].Data(),b);
        corrXF[b] = new TH2D(plotN,plotT,
          N_P_BINS, -1, 1, N_P_BINS, -1, 1);

        plotT = obsT[o] + " #phi_{h} vs. p, for " + cutT + 
          ";p [GeV];#phi_{h}";
        plotN = Form("%s_PhiHvsP_%d",obsN[o].Data(),b);
        PhiHvsP[b] = new TH2D(plotN,plotT,
          N_P_BINS, 0.5, 3*pMaxLo, N_P_BINS, -PI, PI);
        Tools::BinLog(PhiHvsP[b]->GetXaxis());

        plotT = obsT[o] + " #phi_{R} vs. p, for " + cutT + 
          ";p [GeV];#phi_{R}";
        plotN = Form("%s_PhiRvsP_%d",obsN[o].Data(),b);
        PhiRvsP[b] = new TH2D(plotN,plotT,
          N_P_BINS, 0.5, 3*pMaxLo, N_P_BINS, -PI, PI);
        Tools::BinLog(PhiRvsP[b]->GetXaxis());

        plotT = obsT[o] + " #theta vs. p, for " + cutT + 
          ";p [GeV];#theta";
        plotN = Form("%s_thetaVsP_%d",obsN[o].Data(),b);
        thetaVsP[b] = new TH2D(plotN,plotT,
          N_P_BINS, 0.5, 3*pMaxLo, N_P_BINS, 0, PI);
        Tools::BinLog(thetaVsP[b]->GetXaxis());

        plotT = obsT[o] + " q_{T} vs. y, for " + cutT + 
          ";y;q_{T} [GeV]";
        plotN = Form("%s_QtVsY_%d",obsN[o].Data(),b);
        QtVsY[b] = new TH2D(plotN,plotT,
          N_P_BINS, 1e-3, 1,
          N_P_BINS, ptMin, 100*ptMax);
        Tools::BinLog(QtVsY[b]->GetXaxis());
        Tools::BinLog(QtVsY[b]->GetYaxis());

        plotT = obsT[o] + " M_{h} distribution, for " + cutT + ";M_{h} [GeV]";
        plotN = Form("%s_Mh_%d",obsN[o].Data(),b);
        distMh[b] = new TH1D(plotN,plotT, 3*N_P_BINS, 0, 3);
        distMh[b]->SetFillColor(kRed-3);
        distMh[b]->SetLineColor(kRed-3);

        plotT = obsT[o] + " M_{X} distribution, for " + cutT + ";M_{X} [GeV]";
        plotN = Form("%s_Mx_%d",obsN[o].Data(),b);
        distMx[b] = new TH1D(plotN,plotT,
          3*N_P_BINS, conf->bdMmiss[0], conf->bdMmiss[1]);
        distMx[b]->SetFillColor(kRed-5);
        distMx[b]->SetLineColor(kRed-5);

        plotT = obsT[o] + " #theta_{CM} distribution, for " + cutT +
          ";#theta_{CM}";
        plotN = Form("%s_theta_%d",obsN[o].Data(),b);
        distTheta[b] = new TH1D(plotN,plotT, 3*N_P_BINS, 0, 3);
        distTheta[b]->SetFillColor(kCyan+3);
        distTheta[b]->SetLineColor(kCyan+3);
      };
      
      if(o==kEle) {

        plotT = obsT[o] + " y distribution, for " + cutT + ";y";
        plotN = Form("%s_Y_%d",obsN[o].Data(),b);
        distY[b] = new TH1D(plotN,plotT, 3*N_P_BINS, 1e-3, 1);
        distY[b]->SetFillColor(kOrange-7);
        distY[b]->SetLineColor(kOrange-7);
        Tools::BinLog(distY[b]->GetXaxis());

        plotT = obsT[o] + " W distribution, for " + cutT + ";W [GeV]";
        plotN = Form("%s_W_%d",obsN[o].Data(),b);
        distW[b] = new TH1D(plotN,plotT,
          3*N_P_BINS, conf->bdW[0], conf->bdW[1]);
        distW[b]->SetFillColor(kViolet+2);
        distW[b]->SetLineColor(kViolet+2);
      };
    };
  };


  // ----------------------------------------------------



  // prepare for event loop
  Bool_t fillPlots;
  Float_t oP,oPt,oTheta,oEta,oZ,Qt;
  for(b=0; b<NBINS; b++) numEvents[b]=0;


  ///////////////////////////////////////////////
  // EVENT LOOP
  ///////////////////////////////////////////////
  printf("begin loop through %lld events...\n",ev->ENT);
  for(int i=0; i<ev->ENT; i++) {
    //if(i>10000) break; // limiter
    ev->GetEvent(i);

    // fill (x,Q2) plots
    Q2vsXfull->Fill(ev->x,ev->Q2);
    Q2vsQ2pythia->Fill(ev->Q2_pythia,ev->Q2);
    XvsXpythia->Fill(ev->x_pythia,ev->x);
    if(ev->Valid()) {
      Q2vsXcut->Fill(ev->x,ev->Q2);
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
              oZ = UNDEF;
              break;
            case kDih:
              oP = ev->Ph;
              oPt = ev->PhPt;
              oTheta = ev->PhTheta;
              oEta = ev->PhEta;
              oZ = ev->Zpair;
              break;
            case kHadA:
              oP = ev->hadP[qA];
              oPt = ev->hadPt[qA];
              oTheta = ev->hadTheta[qA];
              oEta = ev->hadEta[qA];
              oZ = ev->Z[qA];
              break;
            case kHadB:
              oP = ev->hadP[qB];
              oPt = ev->hadPt[qB];
              oTheta = ev->hadTheta[qB];
              oEta = ev->hadEta[qB];
              oZ = ev->Z[qB];
              break;
          };

          accPolarLoP[o][b]->Fill(oTheta,oP);
          accPolarHiP[o][b]->Fill(oTheta,oP);
          PtVsPzLoP[o][b]->Fill(oP*TMath::Cos(oTheta),oPt);
          PtVsPzHiP[o][b]->Fill(oP*TMath::Cos(oTheta),oPt);
          EtaVsP[o][b]->Fill(oP,oEta);
          EtaVsPt[o][b]->Fill(oPt,oEta);
          EtaVsZ[o][b]->Fill(oZ,oEta);
          PhPerpVsPt[o][b]->Fill(oPt,ev->PhPerp);

          Qt = ev->PhPerp / ev->Zpair;
          QtVsPt[o][b]->Fill(oPt,Qt);
          QtVsEta[o][b]->Fill(oEta,Qt);
          PtVsY[o][b]->Fill(ev->y,oPt);

          if(o==kDih) {
            PhiHvsPhiR[b]->Fill(ev->PhiR,ev->PhiH);
            corrXF[b]->Fill(ev->hadXF[qB],ev->hadXF[qA]);
            PhiHvsP[b]->Fill(oP,ev->PhiH);
            PhiRvsP[b]->Fill(oP,ev->PhiR);
            thetaVsP[b]->Fill(oP,ev->theta);
            QtVsY[b]->Fill(ev->y,Qt);
            distMh[b]->Fill(ev->Mh);
            distMx[b]->Fill(ev->Mmiss);
            distTheta[b]->Fill(ev->theta);
          };
          if(o==kEle) {
            distW[b]->Fill(ev->W);
            distY[b]->Fill(ev->y);
          };
        };

        numEvents[b]++;

      };
    };
  };
  for(b=0; b<NBINS; b++) {
    printf("numEvents in bin %d: %.f\n",b,numEvents[b]);
  };




  // (x,Q2) plot canvases
  TCanvas * Q2vsXfullCanv = xQ2Canv(Q2vsXfull);
  TCanvas * Q2vsXcutCanv = xQ2Canv(Q2vsXcut);


  // polar plot canvases
  TCanvas * accPolarLoPCanv[NOBS][NBINS];
  TCanvas * accPolarHiPCanv[NOBS][NBINS];
  TCanvas * PtVsPzLoPCanv[NOBS][NBINS];
  TCanvas * PtVsPzHiPCanv[NOBS][NBINS];
  Double_t polarMaxLoP;
  Double_t polarMaxHiP;
  for(o=0; o<NOBS; o++) {
    polarMaxLoP = polarMaxHiP = 0;
    for(b=0; b<NBINS; b++) {
      polarMaxLoP = accPolarLoP[o][b]->GetMaximum() > polarMaxLoP ?
                    accPolarLoP[o][b]->GetMaximum() : polarMaxLoP;
      polarMaxHiP = accPolarHiP[o][b]->GetMaximum() > polarMaxHiP ?
                    accPolarHiP[o][b]->GetMaximum() : polarMaxHiP;
    };
    for(b=0; b<NBINS; b++) {
      accPolarLoPCanv[o][b] = 
        PolarCanv(accPolarLoP[o][b],"pol colz",polarMaxLoP);
      accPolarHiPCanv[o][b] =
        PolarCanv(accPolarHiP[o][b],"pol colz",polarMaxHiP);
      PtVsPzLoPCanv[o][b] =
        PolarCanv(PtVsPzLoP[o][b],"colz",polarMaxLoP);
      PtVsPzHiPCanv[o][b] =
        PolarCanv(PtVsPzHiP[o][b],"colz",polarMaxHiP);
    };
  };


  // matrixify plots
  TCanvas * accPolarLoPMatrix[NOBS];
  TCanvas * accPolarHiPMatrix[NOBS];
  TCanvas * PtVsPzLoPMatrix[NOBS];
  TCanvas * PtVsPzHiPMatrix[NOBS];
  TCanvas * EtaVsPMatrix[NOBS];
  TCanvas * EtaVsPtMatrix[NOBS];
  TCanvas * EtaVsZMatrix[NOBS];
  TCanvas * PhPerpVsPtMatrix[NOBS];
  TCanvas * QtVsPtMatrix[NOBS];
  TCanvas * QtVsEtaMatrix[NOBS];
  TCanvas * PtVsYMatrix[NOBS];
  TCanvas * PhiHvsPhiRMatrix;
  TCanvas * corrXFMatrix;
  TCanvas * PhiHvsPMatrix;
  TCanvas * PhiRvsPMatrix;
  TCanvas * thetaVsPMatrix;
  TCanvas * QtVsYMatrix;
  TCanvas * distMhMatrix;
  TCanvas * distMxMatrix;
  TCanvas * distThetaMatrix;
  TCanvas * distWMatrix;
  TCanvas * distYMatrix;
  for(o=0; o<NOBS; o++) {
    accPolarLoPMatrix[o] = MatrixifyCanv(accPolarLoPCanv[o]);
    accPolarHiPMatrix[o] = MatrixifyCanv(accPolarHiPCanv[o]);
    PtVsPzLoPMatrix[o] = MatrixifyCanv(PtVsPzLoPCanv[o]);
    PtVsPzHiPMatrix[o] = MatrixifyCanv(PtVsPzHiPCanv[o]);
    EtaVsPMatrix[o] = MatrixifyDist2(EtaVsP[o],1,0,1);
    EtaVsPtMatrix[o] = MatrixifyDist2(EtaVsPt[o],1,0,1);
    EtaVsZMatrix[o] = MatrixifyDist2(EtaVsZ[o],0,0,1);
    PhPerpVsPtMatrix[o] = MatrixifyDist2(PhPerpVsPt[o],1,1,1);
    QtVsPtMatrix[o] = MatrixifyDist2(QtVsPt[o],1,1,1);
    QtVsEtaMatrix[o] = MatrixifyDist2(QtVsEta[o],0,1,1);
    PtVsYMatrix[o] = MatrixifyDist2(PtVsY[o],1,1,1);
  };
  PhiHvsPhiRMatrix = MatrixifyDist2(PhiHvsPhiR,0,0,1);
  corrXFMatrix = MatrixifyDist2(corrXF,0,0,1);
  PhiHvsPMatrix = MatrixifyDist2(PhiHvsP,1,0,1);
  PhiRvsPMatrix = MatrixifyDist2(PhiRvsP,1,0,1);
  thetaVsPMatrix = MatrixifyDist2(thetaVsP,1,0,1);
  QtVsYMatrix = MatrixifyDist2(QtVsY,1,1,1);
  distMhMatrix = MatrixifyDist1(distMh,0,0);
  distMxMatrix = MatrixifyDist1(distMx,0,0);
  distThetaMatrix = MatrixifyDist1(distTheta,0,0);
  distWMatrix = MatrixifyDist1(distW,0,0);
  distYMatrix = MatrixifyDist1(distY,1,1);


  // write output
  Q2vsXfullCanv->Write();
  Q2vsXcutCanv->Write();

  for(o=0; o<NOBS; o++) {
    outfile->mkdir(obsN[o]);
    outfile->cd(obsN[o]);
    accPolarLoPMatrix[o]->Write();
    accPolarHiPMatrix[o]->Write();
    PtVsPzLoPMatrix[o]->Write();
    PtVsPzHiPMatrix[o]->Write();
    EtaVsPMatrix[o]->Write();
    EtaVsPtMatrix[o]->Write();
    if(o!=kEle) EtaVsZMatrix[o]->Write();
    PhPerpVsPtMatrix[o]->Write();
    QtVsPtMatrix[o]->Write();
    QtVsEtaMatrix[o]->Write();
    PtVsYMatrix[o]->Write();
    if(o==kDih) {
      PhiHvsPhiRMatrix->Write();
      distMhMatrix->Write();
      distMxMatrix->Write();
      distThetaMatrix->Write();
      corrXFMatrix->Write();
      PhiHvsPMatrix->Write();
      PhiRvsPMatrix->Write();
      thetaVsPMatrix->Write();
      QtVsYMatrix->Write();
    };
    if(o==kEle) {
      distWMatrix->Write();
      distYMatrix->Write();
    };
    outfile->cd("/");
  };

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
    if(o!=kEle) { for(b=0; b<NBINS; b++) EtaVsZ[o][b]->Write(); };
    for(b=0; b<NBINS; b++) PhPerpVsPt[o][b]->Write();
    for(b=0; b<NBINS; b++) QtVsPt[o][b]->Write();
    for(b=0; b<NBINS; b++) PtVsY[o][b]->Write();
    if(o==kDih) {
      for(b=0; b<NBINS; b++) PhiHvsPhiR[b]->Write();
      for(b=0; b<NBINS; b++) distMh[b]->Write();
      for(b=0; b<NBINS; b++) distMx[b]->Write();
      for(b=0; b<NBINS; b++) distTheta[b]->Write();
      for(b=0; b<NBINS; b++) corrXF[b]->Write();
      for(b=0; b<NBINS; b++) PhiHvsP[b]->Write();
      for(b=0; b<NBINS; b++) PhiRvsP[b]->Write();
      for(b=0; b<NBINS; b++) thetaVsP[b]->Write();
      for(b=0; b<NBINS; b++) QtVsY[b]->Write();
    };
    if(o==kEle) {
      for(b=0; b<NBINS; b++) distW[b]->Write();
      for(b=0; b<NBINS; b++) distY[b]->Write();
    };
  };
  outfile->cd("/");

  Q2vsQ2pythia->Write();
  XvsXpythia->Write();

  outfile->Close();

};



// take a TH2D and make a polar plot
TCanvas * PolarCanv(TH2D * dist, TString format, Double_t max) {

  // instantiate canvas
  TString canvN = TString(dist->GetName()) + "_Canv";
  TCanvas * canv = new TCanvas(canvN,canvN,800,700);
  canv->SetLogz();
  Double_t radius = dist->GetYaxis()->GetXmax();

  // instantiate frame (because TPad::DrawFrame does not use unique names)
  //canv->DrawFrame(-radius,0,radius,radius,dist->GetTitle());
  TH1F * frame = new TH1F(TString(canvN+"frame"),dist->GetTitle(),
    1000,-radius,radius);
  FormatDist(frame);
  frame->SetXTitle(dist->GetXaxis()->GetTitle());
  frame->SetYTitle(dist->GetYaxis()->GetTitle());
  frame->SetBit(TH1::kNoStats);
  frame->SetMinimum(0);
  frame->SetMaximum(radius);
  frame->GetYaxis()->SetLimits(0,radius);
  frame->Draw(" ");

  // draw polar histogram
  FormatDist(dist);
  dist->SetMinimum(1);
  dist->SetMaximum(max);
  dist->Draw(TString(format+" same"));

  // draw circles of constant momentum
  const Int_t nCircles = 3;
  TEllipse * circle[nCircles];
  Double_t circleRadius;
  for(int c=0; c<nCircles; c++) {
    circleRadius = radius - radius*((float)c)/nCircles;
    circle[c] = new TEllipse(0,0,circleRadius,circleRadius,0,180);
    circle[c]->SetFillStyle(3955); // transparent (SetFillColorAlpha fails)
    circle[c]->SetLineWidth(2);
    circle[c]->Draw();
  };

  // draw spokes of constant eta
  const Int_t etaSpokeMax = 4;
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
  for(int c=0; c<2*etaSpokeMax+1; c++) {
    etaSpoke[c]->SetLineWidth(2);
    etaSpoke[c]->Draw();
  };

  // return canvas
  return canv;
};



// draw bin lines on (x,Q2) plane
TCanvas * xQ2Canv(TH2D * dist) {

  // instantiate canvas
  TString canvN = TString(dist->GetName()) + "_Canv";
  TCanvas * canv = new TCanvas(canvN,canvN,800,700);
  canv->SetLogx();
  canv->SetLogy();
  canv->SetLogz();
  canv->SetGrid(1,1);
  FormatDist(dist);
  dist->Draw("colz");

  // bin lines
  TLine * Q2line[NBINS_Q2];
  TLine * xline[NBINS_x];
  for(b=0; b<NBINS_Q2; b++) {
    Q2line[b] = new TLine(xmin,binBounds_Q2[b],xmax,binBounds_Q2[b]);
    Q2line[b]->SetLineWidth(2);
    Q2line[b]->Draw();
  };
  for(b=0; b<NBINS_x; b++) {
    xline[b] = new TLine(binBounds_x[b],Q2min,binBounds_x[b],Q2max);
    xline[b]->SetLineWidth(2);
    xline[b]->Draw();
  };

  return canv;
};


// return a matrix of (x,Q2) bins, containing plots
TCanvas * MatrixifyCanv(TCanvas ** canvArr) {
  TString canvN = TString(canvArr[0]->GetName()) + "Matrix";
  canvN = canvN.ReplaceAll("_0_","_");
  TCanvas * canv = new TCanvas(canvN,canvN,400*NBINS_x,300*NBINS_Q2);
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

TCanvas * MatrixifyDist1(TH1D ** distArr,Bool_t logx,Bool_t logy) {
  TCanvas * canvases[NBINS];
  TString canvN;
  Double_t max = 0;
  for(b=0; b<NBINS; b++) {
    max = distArr[b]->GetMaximum() > max ? distArr[b]->GetMaximum() : max;
  };
  for(b=0; b<NBINS; b++) {
    canvN = TString(distArr[b]->GetName()) + "_Canv";
    canvases[b] = new TCanvas(canvN,canvN,800,700);
    if(logx) canvases[b]->SetLogx();
    if(logy) canvases[b]->SetLogy();
    canvases[b]->SetGrid(1,1);
    //distArr[b]->SetMaximum(max);
    FormatDist(distArr[b]);
    distArr[b]->Draw();
  };
  return MatrixifyCanv(canvases);
};

TCanvas * MatrixifyDist2(TH2D ** distArr,Bool_t logx,Bool_t logy,Bool_t logz) {
  TCanvas * canvases[NBINS];
  TString canvN;
  Double_t max = 0;
  for(b=0; b<NBINS; b++) {
    max = distArr[b]->GetMaximum() > max ? distArr[b]->GetMaximum() : max;
  };
  for(b=0; b<NBINS; b++) {
    canvN = TString(distArr[b]->GetName()) + "_Canv";
    canvases[b] = new TCanvas(canvN,canvN,800,700);
    if(logx) canvases[b]->SetLogx();
    if(logy) canvases[b]->SetLogy();
    if(logz) canvases[b]->SetLogz();
    canvases[b]->SetGrid(1,1);
    distArr[b]->SetMaximum(max);
    FormatDist(distArr[b]);
    distArr[b]->Draw("colz");
  };
  return MatrixifyCanv(canvases);
};


// set font sizes, etc
void * FormatDist(TH1 * dist) {
  //dist->GetXaxis()->SetTitleSize(0.06);
  //dist->GetYaxis()->SetTitleSize(0.06);
  //dist->GetXaxis()->SetLabelSize(0.06);
  //dist->GetYaxis()->SetLabelSize(0.06);
  //dist->GetXaxis()->SetTitleOffset(1.2);
  //dist->GetYaxis()->SetTitleOffset(0.9);
  dist->SetTitleSize(0.04);
};

