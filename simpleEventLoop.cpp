#include <cstdlib>
#include <iostream>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "Math/Vector4D.h"


// Clas12Tool
#include "reader.h"
#include "bank.h"
#include "particle.h"
#include "clas12reader.h"



#if HIPO_VERSION == 4

int main(int argc, char** argv) {

   // ARGUMENTS
   TString infileN;
   Bool_t dump = false;
   if(argc<=1) {
     printf("USAGE: %s [hipo4 file] [dump text]\n",argv[0]);
     exit(0);
   };
   if(argc>1) infileN = TString(argv[1]);
   if(argc>2) dump = true;


   // set output file name
   TString outfileN = "simpleTree.root";
   printf("outfileN = %s\n",outfileN.Data());
   TFile * outfile = new TFile(outfileN,"RECREATE");


   // define tree
   TTree * tree = new TTree("tree","tree");
   Int_t evnum,pid,helicity;
   Float_t px,py,pz;
   Float_t hadE[2],hadPt[2];
   tree->Branch("evnum",&evnum,"evnum/I");
   tree->Branch("helicity",&helicity,"helicity/I");
   /*
   tree->Branch("pid",&pid,"pid/I");
   tree->Branch("px",&px,"px/F");
   tree->Branch("py",&py,"py/F");
   tree->Branch("pz",&pz,"pz/F");
   */
   tree->Branch("hadE",hadE,"hadE[2]/F"); // [0=pi+, 1=pi-]
   tree->Branch("hadPt",hadPt,"hadPt[2]/F");

   const Float_t PION_MASS = 0.139571;
   const Int_t PION_PID = 211;
   Int_t hadPid[2] = { PION_PID, -PION_PID };

   // HIPO4 reader
   clas12::clas12reader reader(infileN.Data());

   ROOT::Math::PxPyPzMVector pvec;
   ROOT::Math::PxPyPzMVector pionVec[2];
   int h;

   if(dump) {
     gSystem->RedirectOutput("simple.dat","w");
     gSystem->RedirectOutput(0);
   };



   // EVENT LOOP ----------------------------------------------
   printf("begin event loop...\n");
   Int_t evCount=0;
   Float_t Emax[2];
   Int_t lim = (Int_t) 1e6;
   while(reader.next()==true) {
     if(evCount>lim) { fprintf(stderr,"--- stopping loop at %d events\n",lim); break; };

     evnum = reader.runconfig()->getEvent();
     helicity = reader.event()->getHelicity();


     // simple particle loop
     /*
     for(auto & part : reader.getDetParticles()) {
       pid = part->getPid();
       px = part->par()->getPx();
       py = part->par()->getPy();
       pz = part->par()->getPz();


       tree->Fill();

       if(dump) {
         if(abs(pid)==PION_PID) {
           pvec.SetCoordinates(px,py,pz,PION_MASS);
           gSystem->RedirectOutput("simple.dat","a");
           printf("%d %d %.2f %.2f\n", evnum, pid, pvec.E(), pvec.Pt() );
           //printf("%d %d %.2f %.2f %.2f\n", evnum, pid, px, py, pz );
           gSystem->RedirectOutput(0);
         };
       };
     };
     */


     // look for pi+pi-
     for(h=0; h<2; h++) {
       Emax[h]=-1;
       for(auto & part : reader.getByID(hadPid[h])) {
         px = part->par()->getPx();
         py = part->par()->getPy();
         pz = part->par()->getPz();
         pvec.SetCoordinates(px,py,pz,PION_MASS);
         if(pvec.E() > Emax[h]) {
           pionVec[h].SetCoordinates(px,py,pz,PION_MASS);
           Emax[h] = pionVec[h].E();
         };
       };
     };

     if(Emax[0]>0 && Emax[1]>0) {
       for(h=0; h<2; h++) {
         hadE[h] = pionVec[h].E();
         hadPt[h] = pionVec[h].Pt();
       };
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
