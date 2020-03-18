#include <cstdlib>
#include <iostream>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"

// DihBsa
#include "Constants.h"
#include "Tools.h"
#include "EventTree.h"

TString infiles;
Int_t whichPair;

int main(int argc, char** argv) {

   // ARGUMENTS
   infiles = "outroot";
   whichPair = EncodePairType(kPip,kPim);
   if(argc>1) infiles = TString(argv[1]);
   if(argc>2) whichPair = (Int_t)strtof(argv[2],NULL);

   Double_t count=0;
   

   EventTree * ev = new EventTree(infiles,whichPair);
   for(int i=0; i<ev->ENT; i++) {

     ev->GetEvent(i);

     if(ev->Valid() && ev->SpinState()!=0 && 
        ev->Mh>=0 && ev->Mh<=3) {
       count+=1;
     };

   };

   printf("count = %.0f\n",count);
};
