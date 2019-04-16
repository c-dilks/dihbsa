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



int main(int argc, char** argv) {

   // ARGUMENTS
   TString infileN;
   if(argc>1) infileN = TString(argv[1]);
   else {
     printf("USAGE: %s [hipo file]\n",argv[0]);
     exit(0);
   };


   // if true, print more to stdout
   bool debug = 1;

   // set output file name
   TString outfileN;
   outfileN = "out.root";


   // load libs
   gSystem->Load("src/DihBsa.so");
   DIS * disEv = new DIS();
   Trajectory * hadron[2];
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


   // define reader and particle list
   hipo::reader reader;
   reader.open(infileN.Data());
   clas12::particle particleList("REC::Particle",reader);



   // define observable variables
   // -- local index is for the list of observables used here;
   //    this makes looping easier
   // -- global index is the index from src/Constants.h
   // -- example: PartPID(k[o]) gets PID of observable o
   enum obs_enum {oE, oPip, oPim, Nobs}; // local index of observable
   int k[Nobs]; // global particle index of observable
   k[oE] = kE;
   k[oPip] = kPip;
   k[oPim] = kPim;

   clas12::vector4 vecObs[Nobs]; // observable 4-momentum
   Float_t E[Nobs]; // energy 
   Float_t Emax[Nobs]; // maximum energy 
   Int_t Idx[Nobs]; // index of observable in particleList
   Float_t P[Nobs][3]; // momentum [x,y,z]
   enum xyz {dX,dY,dZ};


   Bool_t foundAllObservables;
   Bool_t disCut;


   
   while(reader.next()==true) {
     bench.resume();

     int size = particleList.getSize();

     // search for highest-energy observables
     // -- reset indices and max energy
     for(int o=0; o<Nobs; o++) {
       Emax[o] = 0;
       Idx[o] = -1;
     };
     // -- loop through full particle list, getting max energies
     //    of each observable
     for(int i=0; i<size; i++) {
       for(int o=0; o<Nobs; o++) {
         if( particleList.getPid(i) == PartPID(k[o]) ) {
           particleList.getVector4( i, vecObs[o], PartMass(k[o]) );
           E[o] = vecObs[o].e();
           if(E[o] > Emax[o]) {
             Emax[o] = E[o];
             Idx[o] = i;
           };
         };
       };
     };
     // -- skip this event if not all observables have been found
     foundAllObservables = true;
     for(int o=0; o<Nobs; o++) {
       if(Idx[o]<0) foundAllObservables=false;
     };
     if(!foundAllObservables) continue;
     if(debug) printf("-------\n");
     // -- set energy and momentum for all observables
     for(int o=0; o<Nobs; o++) { 
       E[o] = Emax[o];

       particleList.getVector4( Idx[o], vecObs[o], PartMass(k[o]) ); 

       P[o][dX] = particleList.getPx(Idx[o]);
       P[o][dY] = particleList.getPy(Idx[o]);
       P[o][dZ] = particleList.getPz(Idx[o]);

       //P[o][dX] = vecObs[o].vect().x(); // faster way ??
       //P[o][dY] = vecObs[o].vect().y();
       //P[o][dZ] = vecObs[o].vect().z();
     };



     // compute DIS kinematics
     disEv->SetElectron(
       P[k[oE]][dX],
       P[k[oE]][dY],
       P[k[oE]][dZ]
     );
     disCut = disEv->Analyse();


     // set pi+ and pi- momenta
     hadron[hP]->SetMomentum(
       P[k[oPip]][dX],
       P[k[oPip]][dY],
       P[k[oPip]][dZ]
     );
     hadron[hM]->SetMomentum(
       P[k[oPim]][dX],
       P[k[oPim]][dY],
       P[k[oPim]][dZ]
     );
     for(int h=0; h<2; h++) {
       hadE[h] = hadron[h]->Vec->E();
       hadPt[h] = hadron[h]->Vec->Pt();
       hadEta[h] = hadron[h]->Vec->Eta();
       hadPhi[h] = hadron[h]->Vec->Phi();
     };
       


     tree->Fill();

     bench.pause();
   };
   reader.showBenchmark();
   printf(" time spend on analysis = %8.5f (ms) events = %8d\n",
     bench.getTime()*1e-9,
     bench.getCounter());
   tree->Write("tree");
};
//### END OF GENERATED CODE
