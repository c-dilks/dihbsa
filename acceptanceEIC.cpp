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

// DihBsa
#include "Constants.h"
#include "Tools.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"
#include "EventTree.h"
#include "Config.h"


TString inDir;
Int_t whichPair;
Int_t whichHad[2];
TString hadTitle[2];
TString dihTitle;
int h;


TCanvas * PolarCanv(TH2D * dist);


int main(int argc, char** argv) {

  // ARGUMENTS
  inDir = "outroot";
  whichPair = EncodePairType(kPip,kPim);
  if(argc>1) inDir = TString(argv[1]);
  if(argc>2) whichPair = (Int_t)strtof(argv[2],NULL);

  // get hadron pair from whichPair; note that in the print out, the 
  // order of hadron 0 and 1 is set by Constants::dihHadIdx
  printf("whichPair = 0x%x\n",whichPair);
  DecodePairType(whichPair,whichHad[qA],whichHad[qB]);
  for(h=0; h<2; h++) {
    hadTitle[h] = PairHadTitle(whichHad[qA],whichHad[qB],h);
    printf("hadron %d:  idx=%d  hadron=%s\n",
        h,dihHadIdx(whichHad[qA],whichHad[qB],h),hadTitle[h].Data());
  };
  dihTitle = PairTitle(whichPair);
  printf("dihadron=%s\n",dihTitle.Data());

  EventTree * ev = new EventTree(TString(inDir+"/*.root"),whichPair);
  Config * conf = new Config();


  TFile * outfile = new TFile("eicPlots.root","RECREATE");


  // binning --------------------------------------------
  const Int_t NBINS_Q2 = 6;
  Float_t binBounds_Q2[NBINS_Q2] = {
    0,
    0.1,
    1,
    10,
    100,
    1000
  };
  const Int_t NBINS_x = 4;
  Float_t binBounds_x[NBINS_x] = {
    0,
    0.33,
    0.66,
    1
  };
  //-----------------------------------------------------

  // remap all bin bounds into a 1-dimensional array
  // note: the "last" bin for each variable is an extra; it is used to
  //       integrate over other dimensions (see printout)
  const Int_t NBINS = NBINS_Q2 * NBINS_x;
  Float_t Q2Bin[NBINS][2];
  Float_t xBin[NBINS][2];
  int b;
  for(int bQ2=0; bQ2<NBINS_Q2; bQ2++) {
    for(int bx=0; bx<NBINS_x; bx++) {
      b = bQ2 + bx*NBINS_Q2;
      //b = bQ2 + bx*NBINS_Q2 + bPh*NBINS_x*NBINS_Q2; // for another dim
      for(int bb=0; bb<2; bb++) {
        Q2Bin[b][bb] = (bQ2<NBINS_Q2-1) ? binBounds_Q2[bQ2+bb] :
                                          binBounds_Q2[bb==0?0:NBINS_Q2-1];
        xBin[b][bb] = (bx<NBINS_x-1) ? binBounds_x[bx+bb] :
                                          binBounds_x[bb==0?0:NBINS_x-1];
      };
    };
  };
  for(b=0; b<NBINS; b++) {
    printf("b=%d\tQ2:[%7.2f,%7.2f]\tx[%7.2f,%7.2f]\n",
      b,Q2Bin[b][0],Q2Bin[b][1],
      xBin[b][0],xBin[b][1]
    );
  };



  // define plots ---------------------------------------
  const Int_t NPLOTBINS = 150;
  TH2D * elePolar[NBINS];
  TH2D * dihPolar[NBINS];
  TH2D * hadPolar[NBINS][2];
  TString plotN,plotT,cutT;

  for(b=0; b<NBINS; b++) {

    // cut title
    cutT = Form("%.1f<Q^{2}<%.1f GeV^{2} and %.2f<x<%.2f",
      Q2Bin[b][0],Q2Bin[b][1],
      xBin[b][0],xBin[b][1]
    );

    // -- electron plots
    plotN = Form("elePolar_%d",b);
    plotT = "e^{-} momentum, for " + cutT +
      ";p_{z} [GeV];p_{T} [GeV]";
    elePolar[b] = new TH2D(plotN,plotT,
      NPLOTBINS,-PI,PI,NPLOTBINS,0,conf->EbeamEn);

    // -- dihadron plots
    plotN = Form("dihPolar_%d",b);
    plotT = dihTitle + " momentum, for " + cutT +
      ";P_{h,z} [GeV];P_{h,T} [GeV]";
    dihPolar[b] = new TH2D(plotN,plotT,
        NPLOTBINS,-PI,PI,NPLOTBINS,conf->bdHadP[0],conf->bdHadP[1]);

    // -- hadron plots
    for(h=0; h<2; h++) {
      plotN = Form("had%dPolar_%d",h+1,b);
      plotT = hadTitle[h] + " momentum, for " + cutT + 
        ";p_{z} [GeV];p_{T} [GeV]";
      hadPolar[b][h] = new TH2D(plotN,plotT,
        NPLOTBINS,-PI,PI,NPLOTBINS,conf->bdHadP[0],conf->bdHadP[1]);
    };

  };


  // ----------------------------------------------------




  ///////////////////////////////////////////////
  // EVENT LOOP
  ///////////////////////////////////////////////
  printf("begin loop through %lld events...\n",ev->ENT);
  for(int i=0; i<ev->ENT; i++) {
    ev->GetEvent(i);

    // loop over (Q2,x) bins
    for(b=0; b<NBINS; b++) {

      // override (Q2,x) cuts in EventTree
      ev->cutQ2 = ev->Q2 > Q2Bin[b][0] && ev->Q2 < Q2Bin[b][1];
      ev->cutX = ev->x > xBin[b][0] && ev->x < xBin[b][1];
      ev->cutDIS = ev->cutQ2 && ev->cutX && ev->cutW && ev->cutY;

      // fill polar plots
      if(ev->Valid()) {
        elePolar[b]->Fill(ev->eleTheta,ev->eleP);
        dihPolar[b]->Fill(ev->PhTheta,ev->Ph);
        for(h=0; h<2; h++) 
          hadPolar[b][h]->Fill(ev->hadTheta[h],ev->hadP[h]);
      };
    };
  };


  // define polar plot canvases
  TCanvas * elePolarCanv[NBINS];
  TCanvas * dihPolarCanv[NBINS];
  TCanvas * hadPolarCanv[NBINS][2];
  for(b=0; b<NBINS; b++) {
    elePolarCanv[b] = PolarCanv(elePolar[b]);
    dihPolarCanv[b] = PolarCanv(dihPolar[b]);
    for(h=0; h<2; h++)
      hadPolarCanv[b][h] = PolarCanv(hadPolar[b][h]);
  };


  // draw matrix of polar plots in bins of (Q2,x)
  TString matrixT;
  TCanvas * elePolarMatrix = 
    new TCanvas("elePolarMatrix","elePolarMatrix",1000,1000);
  TCanvas * dihPolarMatrix = 
    new TCanvas("dihPolarMatrix","dihPolarMatrix",1000,1000);
  TCanvas * hadPolarMatrix[2];
  for(h=0; h<2; h++) {
    matrixT = Form("had%dPolarMatrix",h+1);
    hadPolarMatrix[h] = new TCanvas(matrixT,matrixT,1000,1000);
  };
  elePolarMatrix->Divide(NBINS_Q2,NBINS_x);
  dihPolarMatrix->Divide(NBINS_Q2,NBINS_x);
  for(h=0; h<2; h++)
    hadPolarMatrix[h]->Divide(NBINS_Q2,NBINS_x);
  for(b=0; b<NBINS; b++) {
    elePolarMatrix->cd(b+1); elePolarCanv[b]->DrawClonePad();
    dihPolarMatrix->cd(b+1); dihPolarCanv[b]->DrawClonePad();
    for(h=0; h<2; h++) {
      hadPolarMatrix[h]->cd(b+1); hadPolarCanv[b][h]->DrawClonePad();
    };
  };


  // write output
  elePolarMatrix->Write();
  dihPolarMatrix->Write();
  for(h=0; h<2; h++)
    hadPolarMatrix[h]->Write();
  for(b=0; b<NBINS; b++) elePolarCanv[b]->Write();
  for(b=0; b<NBINS; b++) dihPolarCanv[b]->Write();
  for(h=0; h<2; h++) {
    for(b=0; b<NBINS; b++) hadPolarCanv[b][h]->Write();
  };

};



// take a TH2D and make a polar plot
TCanvas * PolarCanv(TH2D * dist) {

  // instantiate canvas
  TString canvN = dist->GetName();
  canvN = canvN.ReplaceAll("Polar","PolarCanv");
  TCanvas * canv = new TCanvas(canvN,canvN,800,700);

  // draw empty distribution, to use its axes
  Double_t radius = dist->GetYaxis()->GetXmax();
  Int_t bgbins = dist->GetYaxis()->GetNbins();
  TString distN = TString(dist->GetName()) + "_bg";
  TString distT = dist->GetTitle();
  TH2D * distBG = new TH2D(distN,distT,bgbins,
    -radius,radius,bgbins/2,0,radius);
  distBG->SetXTitle(dist->GetXaxis()->GetTitle());
  distBG->SetYTitle(dist->GetYaxis()->GetTitle());
  distBG->SetMaximum(dist->GetMaximum());
  distBG->Draw("colz");

  // draw circles of constant momentum
  const Int_t nCircles = 3;
  TEllipse * circle[nCircles+1];
  Double_t circleRadius;
  for(int c=0; c<nCircles; c++) {
    circleRadius = radius - radius*((float)c)/nCircles;
    circle[c] = new TEllipse(0,0,circleRadius,circleRadius,0,180);
    circle[c]->Draw();
  };
  circle[nCircles] = new TEllipse(0,0,1.25,1.25,0,180);
  circle[nCircles]->Draw();

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

  // draw polar plot and return canvas
  dist->Draw("pol colz same");
  return canv;
};
