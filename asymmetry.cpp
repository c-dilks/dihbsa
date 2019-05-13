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
#include "Asymmetry.h"
#include "KinDep.h"

void HadronCompareCanv(TCanvas * canv, TH1F * dist[2], TH2F * corr);

int main(int argc, char** argv) {

   // ARGUMENTS
   TString inDir = "outroot";
   Int_t whichPhiR = 3; // 1:phiRq  2:phiRp_r  3:phiRp
   if(argc>1) inDir = TString(argv[1]);
   if(argc>2) whichPhiR = (Int_t)strtof(argv[2],NULL);

   TFile * outfile = new TFile("spin.root","RECREATE");


   EventTree * ev = new EventTree(TString(inDir+"/*.root"));
   Asymmetry * asym = new Asymmetry(Asymmetry::modSinPhiR, false);
   //Asymmetry * asym = new Asymmetry(Asymmetry::modSinPhiHR, false);


   printf("begin loop through %lld events...\n",ev->ENT);
   for(int i=0; i<ev->ENT; i++) {

     ev->GetEvent(i);

     if(ev->cutDihadron && ev->cutQ2 && ev->cutW && ev->cutY) {
       asym->Mh = ev->Mh;
       asym->x = ev->x;
       asym->z = ev->Zpair;
       asym->eSpin = ev->helicity;
       asym->pSpin = 0;
       asym->PhiH = ev->PhiH;
       asym->PhPerp = ev->PhPerp;

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

       //ev->PrintEvent();
       asym->FillPlots();

     };
   };

   asym->CalculateAsymmetries();

   KinDep * kd = new KinDep(asym);

   asym->WriteObjects(outfile);
   kd->WriteObjects(outfile);
   kd->PrintPNGs(whichPhiR);

   outfile->Close();

};


