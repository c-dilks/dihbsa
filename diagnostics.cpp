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

// DihBsa
#include "Constants.h"
#include "Tools.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"
#include "EventTree.h"

void HadronCompareCanv(TCanvas * canv, TH1D * dist[2], TH2D * corr);
void Hadron2dCanv(TCanvas * canv, TH2D * distA, TH2D * distB);
TString corrTitle(TString var);
TString distTitle(TString var);
TString dist2Title(TString hadron, TString varX,TString varY);

TString inDir;
Int_t whichPair;
Int_t whichHad[2];
TString hadName[2];
TString hadTitle[2];

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
   for(int h=0; h<2; h++) {
     hadName[h] = PairHadName(whichHad[qA],whichHad[qB],h);
     hadTitle[h] = PairHadTitle(whichHad[qA],whichHad[qB],h);
     printf("hadron %d:  idx=%d  name=%s  title=%s\n",
       h,dihHadIdx(whichHad[qA],whichHad[qB],h),hadName[h].Data(),hadTitle[h].Data());
   };

   EventTree * ev = new EventTree(TString(inDir+"/*.root"),whichPair);


   TFile * outfile = new TFile("plots.root","RECREATE");

   const Int_t NBINS = 100;
   Float_t deltaPhi;


   // DIS kinematics
   TH1D * WDist = new TH1D("WDist","W distribution (w/o W cut);W",NBINS,0,6);
   TH1D * XDist = new TH1D("XDist","x distribution;x",NBINS,0,1);
   TH2D * Q2vsW = new TH2D("Q2vsW","Q^{2} vs. W (w/o W cut);W;Q^{2}",
                                   NBINS,0,6,NBINS,0,12);
   TH2D * Q2vsX = new TH2D("Q2vsX","Q^{2} vs. x;x;Q^{2}",NBINS,0,1,NBINS,0,12);
   TH1D * YDist = new TH1D("YDist","y distribution (w/o y cut);y",NBINS,0,1);
   
   // electron kinematics
   TH1D * eleEDist = new TH1D("eleEDist","e^{-} E distribution",NBINS,0,12);
   TH1D * elePtDist = new TH1D("elePtDist","e^{-} p_{T} distribution",NBINS,0,4);
   TH1D * eleEtaDist = new TH1D("eleEtaDist","e^{-} #eta distribution",NBINS,-3,6);
   TH1D * elePhiDist = new TH1D("elePhiDist","e^{-} #phi distribution",NBINS,-PIe,PIe);
   TH2D * eleEtaVsPhi = new TH2D("eleEtavsPhi","e^{-} #eta vs #phi;#phi;#eta",
     NBINS,-PIe,PIe,NBINS,-3,6);
   TH2D * eleEVsPhi = new TH2D("eleEvsPhi","e^{-} E vs #phi;#phi;E",
     NBINS,-PIe,PIe,NBINS,0,12);
   TH2D * elePtVsPhi = new TH2D("elePtvsPhi","e^{-} p_{T} vs #phi;#phi;#p_{T}",
     NBINS,-PIe,PIe,NBINS,0,4);

   // dihadron's hadron kinematic correlations
   TH2D * hadECorr = new TH2D("hadECorr",corrTitle("E"),NBINS,0,10,NBINS,0,10);
   TH2D * hadPCorr = new TH2D("hadPCorr",corrTitle("p"),NBINS,0,10,NBINS,0,10);
   TH2D * hadPtCorr = new TH2D("hadPtCorr",corrTitle("p_{T}"),NBINS,0,4,NBINS,0,4);
   TH2D * hadEtaCorr = new TH2D("hadEtaCorr",corrTitle("#eta"),NBINS,0,5,NBINS,0,5);
   TH2D * hadPhiCorr = new TH2D("hadPhiCorr",corrTitle("#phi"),
                                             NBINS,-PIe,PIe,NBINS,-PIe,PIe);
   TH2D * hadZCorr = new TH2D("hadZCorr",corrTitle("z"),NBINS,0,1,NBINS,0,1);
   
   // dihadron's hadron kinematics
   TH1D * hadEDist[2];
   TH1D * hadPDist[2];
   TH1D * hadPtDist[2];
   TH1D * hadEtaDist[2];
   TH1D * hadPhiDist[2];
   TH1D * hadZDist[2];
   TH2D * hadEtaVsPhi[2];
   TH2D * hadEVsPhi[2];
   TH2D * hadPtVsPhi[2];
   for(int h=0; h<2; h++) {
     hadEDist[h] = new TH1D(TString(hadName[h]+"hadEDist"),distTitle("E"),
       NBINS,0,10);
     hadPDist[h] = new TH1D(TString(hadName[h]+"hadPDist"),distTitle("p"),
       NBINS,0,10);
     hadPtDist[h] = new TH1D(TString(hadName[h]+"hadPtDist"),distTitle("p_{T}"),
       NBINS,0,4);
     hadEtaDist[h] = new TH1D(TString(hadName[h]+"hadEtaDist"),distTitle("#eta"),
       NBINS,0,5);
     hadPhiDist[h] = new TH1D(TString(hadName[h]+"hadPhiDist"),distTitle("#phi"),
       NBINS,-PIe,PIe);
     hadZDist[h] = new TH1D(TString(hadName[h]+"hadZDist"),distTitle("z"),
       NBINS,0,1);

     hadEtaVsPhi[h] = new TH2D(
       TString(hadName[h]+"hadEtaVsPhi"),dist2Title(hadTitle[h],"#phi","#eta"),
       NBINS,-PIe,PIe,NBINS,0,5);
     hadEVsPhi[h] = new TH2D(
       TString(hadName[h]+"hadEVsPhi"),dist2Title(hadTitle[h],"#phi","E"),
       NBINS,-PIe,PIe,NBINS,0,10);
     hadPtVsPhi[h] = new TH2D(
       TString(hadName[h]+"hadPtVsPhi"),dist2Title(hadTitle[h],"#phi","p_{T}"),
       NBINS,-PIe,PIe,NBINS,0,4);
   };


   // dihadron kinematics
   TString plotTitle = "#Delta#phi = #phi(" + hadTitle[qA] + ")" +
                                 " - #phi(" + hadTitle[qB] + 
                                 ") distribution;#Delta#phi";
   TH1D * deltaPhiDist = new TH1D("deltaPhiDist",plotTitle,NBINS,-PIe,PIe);

   TH1D * MhDist = new TH1D("MhDist","M_{h} distribution;M_{h}",2*NBINS,0,3);
   TH1D * ZpairDist = new TH1D("ZpairDist","z_{pair} distribution;z_{pair}",NBINS,0,1);
   TH1D * zetaDist = new TH1D("zetaDist","#zeta distribution;#zeta",NBINS,-1,1);
   TH1D * xFDist = new TH1D("xFDist","x_{F} distribution;x_{F}",NBINS,-2,2);
   TH1D * MmissDist = new TH1D("MmissDist","M_{X} distribution;M_{X}",NBINS,-2,6);
   
   TH1D * PhiHDist = new TH1D("PhiHDist","#phi_{h} distribution;#phi_{h}",
     NBINS,-PIe,PIe);
   TH1D * PhiRDist = new TH1D("PhiRDist","#phi_{R} distribution;#phi_{R}",
     NBINS,-PIe,PIe);
   TH2D * PhiHvsPhiR = new TH2D("PhiHvsPhiR","#phi_{h} vs. #phi_{R};#phi_{R};#phi_{h}",
     NBINS,-PIe,PIe,NBINS,-PIe,PIe);
   TH1D * PhiHRDist = new TH1D("PhiHRDist",
     "#phi_{h}-#phi_{R} distribution;#phi_{h}-#phi_{R}",
     NBINS,-PIe,PIe);

   plotTitle = "P_{h}^{perp}/M_{h} vs. sin(#phi_{h}-#phi_{R});";
   plotTitle += "sin(#phi_{h}-#phi_{R});P_{h}^{perp}/M_{h}";
   TH2D * g1perpWeightVsMod = new TH2D("g1perpWeightVsMod",plotTitle,
     NBINS,-1.1,1.1,
     NBINS,0,6);
   TH2D * PhPerpVsMh = new TH2D("PhPerpVsMh",
     "P_{h}^{perp} vs. M_{h};M_{h};P_{h}^{perp}",
     NBINS,0,3,
     NBINS,0,3);

   // distributions for partial wave analysis
   TH1D * thetaDist = new TH1D("thetaDist","#theta distribution;#theta",NBINS,0,PI);

   TH2D * thetaVsPhiH = new TH2D("thetaVsPhiH",
     "#theta vs #phi_{h};#phi_{h};#theta",
     NBINS,-PIe,PIe,NBINS,0,PIe);
   TH2D * thetaVsPhiR = new TH2D("thetaVsPhiR",
     "#theta vs #phi_{R};#phi_{R};#theta",
     NBINS,-PIe,PIe,NBINS,0,PIe);
   TH2D * thetaVsPhiHR = new TH2D("thetaVsPhiHR",
     "#theta vs #phi_{h}-#phi_{R};#phi_{h}-#phi_{R};#theta",
     NBINS,-PIe,PIe,NBINS,0,PIe);

   TH2D * thetaVsMh = new TH2D("thetaVsMh","#theta vs. M_{h};M_{h};#theta",
     NBINS,0,3,NBINS,0,PIe);
   TH2D * thetaVsZpair = new TH2D("thetaVsZpair","#theta vs. z;z;#theta",
     NBINS,0,1,NBINS,0,PIe);
   TH2D * thetaVsZeta = new TH2D("thetaVsZeta","#theta vs. #zeta;#zeta;#theta",
     NBINS,-1,1,NBINS,0,PIe);
   TH2D * thetaVsX = new TH2D("thetaVsX","#theta vs. x;x;#theta",
     NBINS,0,1,NBINS,0,PIe);
   TH2D * thetaVsPh = new TH2D("thetaVsPh","#theta vs. p;p;#theta",
     NBINS,0,10,NBINS,0,PIe);

   TH2D * thetaVsZ[2];
   TH2D * thetaVsHadP[2];
   for(int h=0; h<2; h++) {
     thetaVsZ[h] = new TH2D(TString("thetaVsZ_"+hadName[h]),
       TString("#theta vs. "+hadTitle[h]+" z;z;#theta"),
       NBINS,0,1,NBINS,0,PIe);
     thetaVsHadP[h] = new TH2D(TString("thetaVsHadP_"+hadName[h]),
       TString("#theta vs. "+hadTitle[h]+" p;p;#theta"),
       NBINS,0,10,NBINS,0,PIe);
   };
   
   TH1D * sinThetaDist = new TH1D("sinThetaDist",
     "sin(#theta) distribution;sin(#theta)",NBINS,-1.1,1.1);
   TH1D * sinThetaCosThetaDist = new TH1D("sinThetaCosThetaDist",
     "sin(#theta)cos(#theta) distribution;sin(#theta)cos(#theta)",NBINS,-1.1,1.1);

   TH1D * cosThetaDist = new TH1D("cosThetaDist",
     "cos(#theta) distribution;cos(#theta)",NBINS,-1.1,1.1);


   // PhiH and PhiR vs. other variables
   TH2D * PhiHvsMh = new TH2D("PhiHvsMh","#phi_{h} vs. M_{h};M_{h};#phi_{h}",
     NBINS,0,3,NBINS,-PIe,PIe);
   TH2D * PhiHvsZ = new TH2D("PhiHvsZ","#phi_{h} vs. z;z;#phi_{h}",
     NBINS,0,1,NBINS,-PIe,PIe);
   TH2D * PhiHvsX = new TH2D("PhiHvsX","#phi_{h} vs. x;x;#phi_{h}",
     NBINS,0,1,NBINS,-PIe,PIe);

   TH2D * PhiRvsMh = new TH2D("PhiRvsMh","#phi_{R} vs. M_{h};M_{h};#phi_{R}",
     NBINS,0,3,NBINS,-PIe,PIe);
   TH2D * PhiRvsZ = new TH2D("PhiRvsZ","#phi_{R} vs. z;z;#phi_{R}",
     NBINS,0,1,NBINS,-PIe,PIe);
   TH2D * PhiRvsX = new TH2D("PhiRvsX","#phi_{R} vs. x;x;#phi_{R}",
     NBINS,0,1,NBINS,-PIe,PIe);

   TH2D * PhiRvsAlpha = new TH2D("PhiRvsAlpha",
     "#phi_{R} vs. #alpha;#alpha;#phi_{R}",
     NBINS,0,1.3,NBINS,-PIe,PIe);
   TH2D * PhiHRvsAlpha = new TH2D("PhiHRvsAlpha",
     "#phi_{h}-#phi_{R} vs. #alpha;#alpha;#phi_{h}-#phi_{R}",
     NBINS,0,1.3,NBINS,-PIe,PIe);
     


   // hadron type matrix
   TH2D * hadTypeMatrix = new TH2D("hadTypeMatrix",
     "Dihadron hadron types matrix;hadron 2;hadron 1",
     nObservables,0,nObservables,nObservables,0,nObservables);
   for(int h1=0; h1<nObservables; h1++) {
     for(int h2=0; h2<nObservables; h2++) {
       hadTypeMatrix->GetYaxis()->SetBinLabel(h1+1,ObsTitle(h1));
       hadTypeMatrix->GetXaxis()->SetBinLabel(h2+1,ObsTitle(h2));
     };
   };


   // multiplicities
   TH1D * partMultiplicity = new TH1D("partMultiplicity",
     "overall particle multiplicities (DIS cuts only)",nParticles,0,nParticles);
   TH1D * obsMultiplicity = new TH1D("obsMultiplicity",
     "dihadrons' particle multiplicities",nObservables,0,nObservables);
   for(int p=0; p<nParticles; p++) 
     partMultiplicity->GetXaxis()->SetBinLabel(p+1,PartTitle(p));
   for(int p=0; p<nObservables; p++) 
     obsMultiplicity->GetXaxis()->SetBinLabel(p+1,ObsTitle(p));
   

   // fiducial phi mask
   TH1D * fiducialPhiMask = new TH1D("fiducialPhiMask",
     "#phi fiducial regions",NBINS,-PIe,PIe);


   // event-level distributions
   TH1D * helicityDist = new TH1D("helicityDist","helicity",5,-2,3);
   TH1D * torusDist = new TH1D("torusDist","torus",5,-2,3);
  
   // diphoton-relevant distributions
   TH1D * diphMdist = new TH1D("diphMdist",
     "M_{#gamma#gamma} distribution;M_{#gamma#gamma}",NBINS,0,1);






   printf("begin loop through %lld events...\n",ev->ENT);
   Int_t hadI[2];
   for(int i=0; i<ev->ENT; i++) {

     ev->GetEvent(i);


     // fill multiplicity plots
     //------------------------
     // fill overall particle multiplicity
     if(ev->cutDIS) {
       for(int p=0; p<nParticles; p++) {
         if(ev->particleCnt[p]>0) partMultiplicity->Fill(p,ev->particleCnt[p]);
       };
     };
     if(ev->cutDIS && ev->cutDihadronKinematics && ev->cutDiph[qA] && ev->cutDiph[qB]) {

       // fill observable multiplicity
       for(int h=0; h<2; h++) {
         hadI[h] = IO(ev->hadIdx[h]); // observable indices
         obsMultiplicity->Fill(hadI[h]);
       };
       
       // fill dihadron matrix
       hadTypeMatrix->Fill(hadI[qB],hadI[qA]);
       // fill transpose elements too (makes symmetric matrix), but don't double-fill
       // diagonal elements
       if(hadI[qA]!=hadI[qB]) hadTypeMatrix->Fill(hadI[qA],hadI[qB]);

     };


     // fill DIS kinematic plots
     // ------------------------
     if(ev->cutDihadron) {

       //ev->PrintEvent();

       if(ev->cutQ2 && ev->cutY) {
         WDist->Fill(ev->W);
         Q2vsW->Fill(ev->W,ev->Q2);
       };

       if(ev->cutQ2 && ev->cutW) {
         YDist->Fill(ev->y);
       };
     };


     // fill dihadron kinematics plots
     // ------------------------------
     if(ev->Valid()) {

       eleEDist->Fill(ev->eleE);
       elePtDist->Fill(ev->elePt);
       eleEtaDist->Fill(ev->eleEta);
       elePhiDist->Fill(ev->elePhi);
       eleEtaVsPhi->Fill(ev->elePhi,ev->eleEta);
       eleEVsPhi->Fill(ev->elePhi,ev->eleE);
       elePtVsPhi->Fill(ev->elePhi,ev->elePt);

       if(Tools::PhiFiducialCut(ev->elePhi)) fiducialPhiMask->Fill(ev->elePhi);


       Q2vsX->Fill(ev->x,ev->Q2);
       XDist->Fill(ev->x);


       hadECorr->Fill(ev->hadE[qB],ev->hadE[qA]);
       hadPCorr->Fill(ev->hadP[qB],ev->hadP[qA]);
       hadPtCorr->Fill(ev->hadPt[qB],ev->hadPt[qA]);
       hadEtaCorr->Fill(ev->hadEta[qB],ev->hadEta[qA]);
       hadPhiCorr->Fill(ev->hadPhi[qB],ev->hadPhi[qA]);
       hadZCorr->Fill(ev->Z[qB],ev->Z[qA]);

       for(int h=0; h<2; h++) {
         hadEDist[h]->Fill(ev->hadE[h]);
         hadPDist[h]->Fill(ev->hadP[h]);
         hadPtDist[h]->Fill(ev->hadPt[h]);
         hadEtaDist[h]->Fill(ev->hadEta[h]);
         hadPhiDist[h]->Fill(ev->hadPhi[h]);
         hadZDist[h]->Fill(ev->Z[h]);

         hadEtaVsPhi[h]->Fill(ev->hadPhi[h],ev->hadEta[h]);
         hadEVsPhi[h]->Fill(ev->hadPhi[h],ev->hadE[h]);
         hadPtVsPhi[h]->Fill(ev->hadPhi[h],ev->hadPt[h]);
       };

       deltaPhi = Tools::AdjAngle(ev->hadPhi[qA] - ev->hadPhi[qB]);
       deltaPhiDist->Fill(deltaPhi);

       MhDist->Fill(ev->Mh);
       ZpairDist->Fill(ev->Zpair);
       zetaDist->Fill(ev->zeta);
       xFDist->Fill(ev->xF);
       MmissDist->Fill(ev->Mmiss);

       PhiHDist->Fill(ev->PhiH);
       PhiRDist->Fill(ev->PhiR);
       PhiHvsPhiR->Fill(ev->PhiR,ev->PhiH);

       PhiHRDist->Fill(ev->PhiHR);
       g1perpWeightVsMod->Fill(TMath::Sin(ev->PhiHR),ev->PhPerp/ev->Mh);
       PhPerpVsMh->Fill(ev->Mh,ev->PhPerp);

       thetaDist->Fill(ev->theta);
       sinThetaDist->Fill(TMath::Sin(ev->theta));
       sinThetaCosThetaDist->Fill(TMath::Sin(ev->theta)*TMath::Cos(ev->theta));
       cosThetaDist->Fill(TMath::Cos(ev->theta));

       thetaVsPhiH->Fill(ev->PhiH,ev->theta);
       thetaVsPhiR->Fill(ev->PhiR,ev->theta);
       thetaVsPhiHR->Fill(ev->PhiHR,ev->theta);

       thetaVsMh->Fill(ev->Mh,ev->theta);
       thetaVsX->Fill(ev->x,ev->theta);
       thetaVsZpair->Fill(ev->Zpair,ev->theta);
       thetaVsZeta->Fill(ev->zeta,ev->theta);
       thetaVsPh->Fill(ev->Ph,ev->theta);
       for(int h=0; h<2; h++) {
         thetaVsZ[h]->Fill(ev->Z[h],ev->theta);
         thetaVsHadP[h]->Fill(ev->hadP[h],ev->theta);
       };

       PhiHvsMh->Fill(ev->Mh,ev->PhiH);
       PhiHvsX->Fill(ev->x,ev->PhiH);
       PhiHvsZ->Fill(ev->Zpair,ev->PhiH);

       PhiRvsMh->Fill(ev->Mh,ev->PhiR);
       PhiRvsX->Fill(ev->x,ev->PhiR);
       PhiRvsZ->Fill(ev->Zpair,ev->PhiR);
       PhiRvsAlpha->Fill(ev->alpha,ev->PhiR);
       PhiHRvsAlpha->Fill(ev->alpha,ev->PhiHR);

       helicityDist->Fill(ev->helicity);
       torusDist->Fill(ev->torus);

       for(int dp=0; dp<ev->diphCnt; dp++) diphMdist->Fill(ev->diphM[dp]);

     };


   }; // eo event loop


   WDist->Write();
   XDist->Write();
   Q2vsW->Write();
   Q2vsX->Write();
   YDist->Write();

   eleEDist->Write();
   elePtDist->Write();
   eleEtaDist->Write();
   elePhiDist->Write();
   eleEtaVsPhi->Write();
   eleEVsPhi->Write();
   elePtVsPhi->Write();
   fiducialPhiMask->Write();



   partMultiplicity->Write();
   obsMultiplicity->Write();
   hadTypeMatrix->Write();


   TCanvas * hadECanv = new TCanvas("hadECanv","hadECanv",1000,800);
   TCanvas * hadPCanv = new TCanvas("hadPCanv","hadPCanv",1000,800);
   TCanvas * hadPtCanv = new TCanvas("hadPtCanv","hadPtCanv",1000,800);
   TCanvas * hadEtaCanv = new TCanvas("hadEtaCanv","hadEtaCanv",1000,800);
   TCanvas * hadPhiCanv = new TCanvas("hadPhiCanv","hadPhiCanv",1000,800);
   TCanvas * hadZCanv = new TCanvas("hadZCanv","hadZCanv",1000,800);
   TCanvas * hadEtaVsPhiCanv = new TCanvas("hadEtaVsPhiCanv","hadEtaVsPhiCanv",1000,800);
   TCanvas * hadEVsPhiCanv = new TCanvas("hadEVsPhiCanv","hadEVsPhiCanv",1000,800);
   TCanvas * hadPtVsPhiCanv = new TCanvas("hadPtVsPhiCanv","hadPtVsPhiCanv",1000,800);

   HadronCompareCanv(hadECanv, hadEDist, hadECorr);
   HadronCompareCanv(hadPCanv, hadPDist, hadPCorr);
   HadronCompareCanv(hadPtCanv, hadPtDist, hadPtCorr);
   HadronCompareCanv(hadEtaCanv, hadEtaDist, hadEtaCorr);
   HadronCompareCanv(hadPhiCanv, hadPhiDist, hadPhiCorr);
   HadronCompareCanv(hadZCanv, hadZDist, hadZCorr);
   Hadron2dCanv(hadEtaVsPhiCanv, hadEtaVsPhi[qA], hadEtaVsPhi[qB]);
   Hadron2dCanv(hadEVsPhiCanv, hadEVsPhi[qA], hadEVsPhi[qB]);
   Hadron2dCanv(hadPtVsPhiCanv, hadPtVsPhi[qA], hadPtVsPhi[qB]);

   hadECanv->Write();
   hadPCanv->Write();
   hadPtCanv->Write();
   hadEtaCanv->Write();
   hadPhiCanv->Write();
   hadZCanv->Write();
   hadEtaVsPhiCanv->Write();
   hadEVsPhiCanv->Write();
   hadPtVsPhiCanv->Write();


   deltaPhiDist->Write();

   MhDist->Write();
   ZpairDist->Write();
   zetaDist->Write();
   xFDist->Write();
   MmissDist->Write();

   PhiHDist->Write();
   PhiRDist->Write();
   PhiHvsPhiR->Write();

   PhiHRDist->Write();
   g1perpWeightVsMod->Write();
   PhPerpVsMh->Write(); 
  
   thetaDist->Write();
   sinThetaDist->Write();
   cosThetaDist->Write();
   sinThetaCosThetaDist->Write();

   thetaVsPhiH->Write();
   thetaVsPhiR->Write();
   thetaVsPhiHR->Write();

   thetaVsMh->Write();
   thetaVsX->Write();
   thetaVsZpair->Write();
   thetaVsZeta->Write();
   for(int h=0; h<2; h++) thetaVsZ[h]->Write();
   thetaVsPh->Write();
   for(int h=0; h<2; h++) thetaVsHadP[h]->Write();

   PhiHvsMh->Write();
   PhiHvsX->Write();
   PhiHvsZ->Write();

   PhiRvsMh->Write();
   PhiRvsX->Write();
   PhiRvsZ->Write();
   PhiRvsAlpha->Write();
   PhiHRvsAlpha->Write();

   torusDist->Write();
   helicityDist->Write();

   diphMdist->Write();

   outfile->Close();

};



// make title for hadron correlation plots
TString corrTitle(TString var) {
  TString varX = var + "(" + hadTitle[qB] + ")";
  TString varY = var + "(" + hadTitle[qA] + ")";
  TString ret = varY + " vs. " + varX + ";" + varX + ";" + varY;
  return ret;
};


TString distTitle(TString var) {
  TString col[2]; 
  if(whichHad[qA]!=whichHad[qB]) {
    for(int cc=0; cc<2; cc++) {
      col[cc] = PartColorName(dihHadIdx(whichHad[qA],whichHad[qB],cc)) + 
                ":" + hadTitle[cc];
    };
  } else {
    col[qA] = "solid:" + hadTitle[qA];
    col[qB] = "dashed:" + hadTitle[qB]; 
  };
  TString ret = var + " distribution (" + col[qA] + " " + col[qB] + ");" + var;
  return ret;
};


TString dist2Title(TString hadron, TString varX,TString varY) {
  return hadron + " " + varY + " vs. " + varX + ";" + varX + ";" + varY;
};



// make canvas for hadron correlation plots
void HadronCompareCanv(TCanvas * canv, TH1D * dist[2], TH2D * corr) {

  for(int h=0; h<2; h++) {
    dist[h]->SetLineColor(PartColor(dihHadIdx(whichHad[qA],whichHad[qB],h)));
    dist[h]->SetLineWidth(2);
  };
  if(whichHad[qA]==whichHad[qB]) dist[qB]->SetLineStyle(kDashed);

  canv->Divide(2,1);

  canv->cd(1);
  Int_t f = dist[qA]->GetMaximum() > dist[qB]->GetMaximum() ? qA:qB;
  dist[f]->Draw();
  dist[(f+1)%2]->Draw("SAME");

  canv->cd(2);
  corr->Draw("colz");
  canv->GetPad(2)->SetGrid(1,1);
};

// make canvas for 2d plots for each hadron 
void Hadron2dCanv(TCanvas * canv, TH2D * distA, TH2D * distB) {
  canv->Divide(2,1);
  canv->cd(1);
  distA->Draw("colz");
  canv->cd(2);
  distB->Draw("colz");
};

