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
     tree->Branch(TString(partName[h]+"Px"),&Px[h],TString(partName[h]+"Px/F"));
     tree->Branch(TString(partName[h]+"Py"),&Py[h],TString(partName[h]+"Py/F"));
     tree->Branch(TString(partName[h]+"Pz"),&Pz[h],TString(partName[h]+"Pz/F"));
     tree->Branch(TString(partName[h]+"E"),&En[h],TString(partName[h]+"E/F"));
     tree->Branch(TString(partName[h]+"Pt"),&Pt[h],TString(partName[h]+"Pt/F"));
   };


   TLorentzVector pv; // 4-momentum
   Bool_t found[N];


   // define text file
   TString dumpFile = "simple.dat";
   if(dump) {
     gSystem->RedirectOutput(dumpFile,"w");
     gSystem->RedirectOutput(0);
   };



   // define clas12reader instance
   clas12::clas12reader reader(infileN.Data());

   // define extra stuff to access MC banks
   hipo::dictionary factory; 
   reader.getReader().readDictionary(factory);
   hipo::event readerEvent;
   hipo::bank mcParticle(factory.getSchema("MC::Particle"));

   ///////////////////////
   Bool_t useMC = 0; // if true, use MC bank instead of REC::Particle bank
   ///////////////////////


   // EVENT LOOP ----------------------------------------------
   printf("begin event loop...\n");
   Int_t nTotal=0;
   Int_t nFound=0;
   while(reader.next()==true) {
     //if(nTotal>1e5) { fprintf(stderr,"--- stop loop at %d events\n",nTotal); break; };

     for(h=0; h<N; h++) {
       En[h] = -1;
       found[h] = false;
     };


     evnum = reader.runconfig()->getEvent();
     helicity = reader.event()->getHelicity();

     // REC::Particle
     if(!useMC) {
       for(h=0; h<N; h++) {
         for(auto & part : reader.getByID(partPid[h])) {
           pv.SetXYZM(
             part->par()->getPx(),
             part->par()->getPy(),
             part->par()->getPz(),
             partMass[h]
           );

           // print each particle
           printf("%d %d %f %f %f\n",evnum,partPid[h],
             part->par()->getPx(),
             part->par()->getPy(),
             part->par()->getPz()
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
     }

     // MC::Particle
     else if(useMC) {
       reader.getReader().read(readerEvent);
       readerEvent.getStructure(mcParticle);
       for(int rr=0; rr<mcParticle.getRows(); rr++) {
         for(h=0; h<N; h++) {
           if(mcParticle.getInt("pid",rr) == partPid[h]) {
             pv.SetXYZM(
               mcParticle.getFloat("px",rr),
               mcParticle.getFloat("py",rr),
               mcParticle.getFloat("pz",rr),
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
       };
     };


     if(found[kE] && found[kP] && found[kM]) {

       tree->Fill();

       if(dump) {
         gSystem->RedirectOutput(dumpFile,"a");
         printf("%d",evnum);
         for(int j : hadrons) {
           printf(" %.2f %.2f %.2f",Px[j],Py[j],Pz[j]);
           printf(" %.2f %.2f",En[j],Pt[j]);
         };
         printf("\n");
         gSystem->RedirectOutput(0);
       };
       nFound++;
     };

     nTotal++;

   };
   // END EVENT LOOP ------------------------------------------


   // write output tree and close output file
   printf("writing tree...\n");
   tree->Write("tree");
   printf("tree written\n");
   outfile->Close();
   printf("\n%s written\n\n",outfileN.Data());

   printf("%d / %d events had a e,pi+,pi-\n",nFound,nTotal);
};
