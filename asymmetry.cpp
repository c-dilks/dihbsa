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

void HadronCompareCanv(TCanvas * canv, TH1F * dist[2], TH2F * corr);

int main(int argc, char** argv) {

   // ARGUMENTS
   TString inDir = "outroot";
   Bool_t useBreit = false;
   if(argc>1) inDir = TString(argv[1]);
   if(argc>2) useBreit = (Bool_t)strtof(argv[2],NULL);

   EventTree * ev = new EventTree(TString(inDir+"/out*.root"));
   TString frame = useBreit ? " --- Breit frame":" --- Lab frame ";
   printf("\n%s\n",frame.Data());

   Asymmetry * asym = new Asymmetry(Asymmetry::modSinPhiR, false);
   


   TFile * outfile = new TFile("spin.root","RECREATE");



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

       asym->PhiR = ev->PhiRq; ///////////////

       asym->FillPlots();

     };
   };

   asym->CalculateAsymmetries();

   asym->Write(outfile);


   outfile->Close();

};


