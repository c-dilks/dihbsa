#include <cstdlib>
#include <iostream>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"

// Clas12Tool
#include "reader.h"
#include "node.h"
#include "bank.h"
#include "particle.h"
#include "detector.h"

// DihBsa
#include "Constants.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"



int main(int argc, char** argv) {

   // ARGUMENTS
   TString infileN;
   if(argc>1) infileN = TString(argv[1]);
   else {
     printf("USAGE: %s [hipo file]\n",argv[0]);
     exit(0);
   };


   // if true, print more to stdout
   bool debug = 0;

   // set output file name
   TString outfileN;
   outfileN = "out.root";


   // load libs
   gSystem->Load("src/DihBsa.so");
   DIS * disEv = new DIS();
   Trajectory * hadron[2];
   Dihadron * dih = new Dihadron();

   enum plus_minus { hP, hM };
   hadron[hP] = new Trajectory(kPip);
   hadron[hM] = new Trajectory(kPim);

   hipo::benchmark bench;


   // define tree
   TFile * outfile = new TFile(outfileN,"RECREATE");
   TTree * tree = new TTree();

   // DIS kinematics branches
   tree->Branch("W",&(disEv->W),"W/F");
   tree->Branch("Q2",&(disEv->Q2),"Q2/F");
   tree->Branch("Nu",&(disEv->Nu),"Nu/F");
   tree->Branch("X",&(disEv->X),"X/F");

   // hadron branches
   Float_t hadE[2]; // [enum plus_minus (0=+, 1=-)]
   Float_t hadPt[2];
   Float_t hadEta[2];
   Float_t hadPhi[2];
   tree->Branch("E",hadE,"E[2]/F");
   tree->Branch("Pt",hadPt,"Pt[2]/F");
   tree->Branch("Eta",hadEta,"Eta[2]/F");
   tree->Branch("Phi",hadPhi,"Phi[2]/F");

   // dihadron branches
   Float_t Mh;
   tree->Branch("Mh",&Mh,"Mh/F");


   // define reader and particle list
   hipo::reader reader;
   reader.open(infileN.Data());
   clas12::particle particleList("REC::Particle",reader);



   // define observable variables
   Float_t E[nParticles]; // energy 
   Float_t Emax[nParticles]; // maximum energy 
   Int_t Idx[nParticles]; // index of observable in particleList
   Float_t P[nParticles][3]; // momentum [x,y,z]
   enum xyz {dX,dY,dZ};
   clas12::vector4 vecObs; // observable 4-momentum
   Int_t oCur,pidCur;


   Bool_t foundAllObservables;
   Bool_t disCut;

   Int_t evCount = 0;

   printf("begin event loop...");
   while(reader.next()==true) {
     bench.resume();


     // search for highest-energy observables
     // -- reset indices and max energy
     for(int o=0; o<nParticles; o++) {
       Emax[o] = 0;
       Idx[o] = -1;
     };
     // -- if it's an electron or pion, check if it has
     //    higher energy than any previous one found
     for(int i=0; i<particleList.getSize(); i++) {

       pidCur = particleList.getPid(i);
       if      (pidCur == PartPID(kE)   )  oCur=kE;
       else if (pidCur == PartPID(kPip) )  oCur=kPip;
       else if (pidCur == PartPID(kPim) )  oCur=kPim;
       else oCur=-10000;

       if(oCur>-10000) {
         particleList.getVector4(i,vecObs,PartMass(oCur));
         E[oCur] = vecObs.e();
         if(E[oCur] > Emax[oCur]) {
           Emax[oCur] = E[oCur];
           Idx[oCur] = i;
         };
       };
     };
     // -- check if we have an electron and pi+/pi- dihadron in this
     //    event; skip if we don't
     foundAllObservables = Idx[kE]>=0 && Idx[kPip]>=0 && Idx[kPim]>=0;
     if(!foundAllObservables) continue;
     if(debug) printf(">>> BEGIN DIHADRON EVENT\n");
     // -- set energy and momentum for all observables
     for(int o=0; o<nParticles; o++) { 
       if(Idx[o]>=0) {
         E[o] = Emax[o];

         particleList.getVector4(Idx[o],vecObs,PartMass(o)); 

         P[o][dX] = particleList.getPx(Idx[o]);
         P[o][dY] = particleList.getPy(Idx[o]);
         P[o][dZ] = particleList.getPz(Idx[o]);
       };
     };


     // compute DIS kinematics
     disEv->SetElectron(
       P[kE][dX],
       P[kE][dY],
       P[kE][dZ]
     );
     disEv->Analyse();

     //disCut = disEv->W > 2.0  &&  disEv->Q2 > 1.0;
     disCut = disEv->W > 2.0;

     if(!disCut) continue;


     // set pi+ and pi- momenta
     hadron[hP]->SetMomentum(
       P[kPip][dX],
       P[kPip][dY],
       P[kPip][dZ]
     );
     hadron[hM]->SetMomentum(
       P[kPim][dX],
       P[kPim][dY],
       P[kPim][dZ]
     );
     for(int h=0; h<2; h++) {
       hadE[h] = hadron[h]->Vec->E();
       hadPt[h] = hadron[h]->Vec->Pt();
       hadEta[h] = hadron[h]->Vec->Eta();
       hadPhi[h] = hadron[h]->Vec->Phi();
       if(debug) {
         printf("[+] %s 4-momentum:\n",(hadron[h]->Title()).Data());
         hadron[h]->Vec->Print();
       };
     };


     // set dihadron momenta
     dih->SetHadrons(hadron[hP],hadron[hM]);
     Mh = dih->Mh();

       


     tree->Fill();
     evCount++;
     /*
     if(evCount%100==0) 
       printf("[---] %d dihadron events found\n",evCount);
       */

     bench.pause();
   };
   reader.showBenchmark();
   printf(" time spend on analysis = %8.5f (ms) events = %8d\n",
     bench.getTime()*1e-9,
     bench.getCounter());
   tree->Write("tree");
   printf("\nDIHADRON EVENT COUNT: %d\n\n",evCount);
};
//### END OF GENERATED CODE
