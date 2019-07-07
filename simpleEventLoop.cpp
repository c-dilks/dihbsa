#include <cstdlib>
#include <iostream>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLorentzVector.h"


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


   TString outfileN = "simpleTree.root";
   printf("outfileN = %s\n",outfileN.Data());
   TFile * outfile = new TFile(outfileN,"RECREATE");


   // define observables
   enum part_enum { kE, kP, kM, N };
   TString partName[N]; partName[kE] = "ele"; partName[kP] = "pip"; partName[kM] = "pim";
   Int_t partPid[N];    partPid[kE] = 11;     partPid[kP] = 211;    partPid[kM] = -211;
   Float_t partMass[N];
   partMass[kE] = 0.000511;
   partMass[kP] = 0.139571;
   partMass[kM] = 0.139571;
   int h;
   int hadrons[2] = {kP,kM};


   // define tree
   TTree * tree = new TTree("tree","tree");
   Int_t evnum,helicity;
   Float_t En[N],Pt[N],Px[N],Py[N],Pz[N];
   tree->Branch("evnum",&evnum,"evnum/I");
   tree->Branch("helicity",&helicity,"helicity/I");
   for(h=0; h<N; h++) {
     tree->Branch(TString(partName[h]+"E"),&En[h],TString(partName[h]+"E/F"));
     tree->Branch(TString(partName[h]+"Pt"),&Pt[h],TString(partName[h]+"Pt/F"));
     tree->Branch(TString(partName[h]+"Px"),&Px[h],TString(partName[h]+"Px/F"));
     tree->Branch(TString(partName[h]+"Py"),&Py[h],TString(partName[h]+"Py/F"));
     tree->Branch(TString(partName[h]+"Pz"),&Pz[h],TString(partName[h]+"Pz/F"));
   };


   TLorentzVector pv; // 4-momentum
   Bool_t found[N];


   // define text file
   TString dumpFile = "simple.dat";
   if(dump) {
     gSystem->RedirectOutput(dumpFile,"w");
     gSystem->RedirectOutput(0);
   };



   // EVENT LOOP ----------------------------------------------
   clas12::clas12reader reader(infileN.Data());
   printf("begin event loop...\n");
   Int_t evCount=0;
   Int_t lim = (Int_t) 1e6;
   while(reader.next()==true) {
     if(evCount>lim) { fprintf(stderr,"--- stopping loop at %d events\n",lim); break; };

     for(h=0; h<N; h++) {
       En[h] = -1;
       found[h] = false;
     };


     evnum = reader.runconfig()->getEvent();
     helicity = reader.event()->getHelicity();

     for(h=0; h<N; h++) {
       for(auto & part : reader.getByID(partPid[h])) {
         pv.SetXYZM(
           part->par()->getPx(),
           part->par()->getPy(),
           part->par()->getPz(),
           partMass[h]
         );
         if(pv.E() > En[h]) {
           En[h] = pv.E();
           Pt[h] = pv.Pt();
           Px[h] = pv.Px();
           Py[h] = pv.Py();
           Pz[h] = pv.Pz();
           found[h] = true;
         };
       };
     };


     if(found[kE] && found[kP] && found[kM]) {

       tree->Fill();

       if(dump) {
         gSystem->RedirectOutput(dumpFile,"a");
         printf("%d",evnum);
         for(int j : hadrons) {
           //printf(" %.2f %.2f",En[j],Pt[j]);
           printf(" %.2f %.2f %.2f",Px[j],Py[j],Pz[j]);
         };
         printf("\n");
         gSystem->RedirectOutput(0);
       };
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
