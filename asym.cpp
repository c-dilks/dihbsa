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

// DihBsa
#include "Constants.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"
#include "EventTree.h"
#include "Binning.h"
#include "Asymmetry.h"


int main(int argc, char** argv) {

   // ARGUMENTS
   TString inDir = "outroot";
   Int_t whichModulation = 0; // see src/Asymmetry.h
   Int_t whichPhiR = 3; // 1:phiRq  2:phiRp_r  3:phiRp // Alessandro prefers 3:phiRp
   if(argc>1) inDir = TString(argv[1]);
   if(argc>2) whichModulation = (Int_t)strtof(argv[2],NULL);
   if(argc>3) whichPhiR = (Int_t)strtof(argv[2],NULL);

   Binning * B = new Binning();

   TFile * outfile = new TFile("spin.root","RECREATE");


   EventTree * ev = new EventTree(TString(inDir+"/*.root"));
   Asymmetry * asym = new Asymmetry(B, whichModulation, 1, Binning::vM, 0);
   if(!(asym->success)) return 0;

   //outfile->Close(); return 0;


   printf("begin loop through %lld events...\n",ev->ENT);
   for(int i=0; i<ev->ENT; i++) {

     ev->GetEvent(i);

     if(ev->cutDihadron && ev->cutQ2 && ev->cutW && ev->cutY) {

       switch(whichPhiR) {
         case 1:
           asym->PhiR = ev->PhiRq;
           break;
         case 2:
           asym->PhiR = ev->PhiRp_r;
           break;
         case 3:
           asym->PhiR = ev->PhiRp;
           break;
         default:
           fprintf(stderr,"ERROR: invalid whichPhiR\n");
           return 0;
       };

       if(B->GetBin(Binning::vM, ev->Mh) == 0) {

         //ev->PrintEvent();
         
         asym->Mh = ev->Mh;
         asym->x = ev->x;
         asym->z = ev->Zpair;
         asym->eSpin = ev->helicity;
         asym->pSpin = 0;
         asym->PhiH = ev->PhiH;
         asym->PhPerp = ev->PhPerp;
         asym->FillPlots();
       };

     };
   };

   asym->CalculateAsymmetries();


   asym->ivDist1->Write();
   asym->modDist->Write();
   for(Int_t m=0; m<Asymmetry::nModBins; m++) asym->modBinDist[m]->Write();
   asym->IVvsModDist->Write();
   for(Int_t s=0; s<nSpin; s++) asym->aziDist[s]->Write();
   asym->asymGr->Write();

   outfile->Close();

   printf("end %s\n",argv[0]);

};


