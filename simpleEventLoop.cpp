#include <cstdlib>
#include <iostream>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TObjArray.h"


// Clas12Tool
#include "reader.h"
#include "bank.h"
#include "particle.h"
#include "clas12reader.h"

// DihBsa
#include "Constants.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"
#include "Diphoton.h"


#if HIPO_VERSION == 4

int main(int argc, char** argv) {

   // ARGUMENTS
   TString infileN;
   if(argc<=1) {
     printf("USAGE: %s [hipo4 file]\n",argv[0]);
     exit(0);
   };
   if(argc>1) infileN = TString(argv[1]);


   // set output file name
   TString outfileN = "simpleTree.root";
   printf("outfileN = %s\n",outfileN.Data());
   TFile * outfile = new TFile(outfileN,"RECREATE");


   // define tree
   TTree * tree = new TTree();
   Int_t evnum,pid,helicity;
   Float_t px,py,pz;
   tree->Branch("evnum",&evnum,"evnum/I");
   tree->Branch("helicity",&helicity,"helicity/I");
   tree->Branch("pid",&pid,"pid/I");
   tree->Branch("px",&px,"px/F");
   tree->Branch("py",&py,"py/F");
   tree->Branch("pz",&pz,"pz/F");


   // HIPO4 reader
   clas12::clas12reader reader(infileN.Data());



   // EVENT LOOP ----------------------------------------------
   printf("begin event loop...\n");
   Int_t evCount=0;
   Int_t lim = (Int_t) 1e6;
   while(reader.next()==true) {
     if(evCount>lim) { fprintf(stderr,"--- stopping loop at %d events\n",lim); break; };

     evnum = reader.runconfig()->getEvent();
     helicity = reader.event()->getHelicity();


     // particle loop
     for(auto & part : reader.getDetParticles()) {
       pid = part->getPid();
       px = part->par()->getPx();
       py = part->par()->getPy();
       pz = part->par()->getPz();

       tree->Fill();

     };

     evCount++;

   };
   // END EVENT LOOP ------------------------------------------


   // write output tree and close output file
   printf("writing tree...\n");
   tree->Write("tree");
   printf("tree written\n");
   outfile->Close();
   printf("\n%s written\n\n",outfileN.Data());
};


#else
int main(int argc, char** argv) {
  fprintf(stderr,"ERROR: must use HIPO_VERSION == 4\n");
  exit(0);
};
#endif