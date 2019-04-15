#include <cstdlib>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "reader.h"
#include "node.h"
#include "bank.h"
#include "particle.h"
#include "detector.h"


// pdg particle ids
static const int PID_E = 11; // electron
static const int PID_P = 2212; // proton
static const int PID_N = 2112; // neutron
static const int PID_PIP = 211; // pi+
static const int PID_PIM = -211; // pi-
static const int PID_PI0 = 111; // pi0
static const int PID_KP = 321; // K+
static const int PID_KM = -321; // K-
static const int PID_PHOTON = 22; // photon

// pdg particle masses [GeV/c^2]
static const double MASS_E = 0.000511;
static const double MASS_P = 0.938272;
static const double MASS_N = 0.939565;
static const double MASS_PIP = 0.139571;
static const double MASS_PIM = 0.139571;
static const double MASS_PI0 = 0.134977;
static const double MASS_KP = 0.493677;
static const double MASS_KM = 0.493677;
static const double MASS_PHOTON = 0.0;

// beam energy -- where can this be read from the HIPO file ??
static const double BEAM_EN = 10.6; // [GeV]


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


   // define tree
   TFile * outfile = new TFile(outfileN,"RECREATE");
   TTree * tree = new TTree();

   Float_t W,Q2,Nu,X;
   tree->Branch("W",&W,"W/F");
   tree->Branch("Q2",&Q2,"Q2/F");
   tree->Branch("Nu",&Nu,"Nu/F");
   tree->Branch("X",&X,"X/F");


   // define reader and particle list
   hipo::reader reader;
   reader.open(infileN.Data());
   clas12::particle particleList("REC::Particle",reader);


   clas12::vector4 vecBeam( 0.0, 0.0, BEAM_EN, BEAM_EN );
   clas12::vector4 vecTarget( 0.0, 0.0, 0.0, MASS_P );
   clas12::vector4 vecW;
   clas12::vector4 vecQ;

   hipo::benchmark b;

   // define observable variables
   enum obs_enum { kE, kP, kM, Nobs};
   TString Name[Nobs];
   Name[kE] = "e-";
   Name[kP] = "pi+";
   Name[kM] = "pi-";
   Int_t PID[Nobs];
   PID[kE] = PID_E;
   PID[kP] = PID_PIP;
   PID[kM] = PID_PIM;
   Double_t Mass[Nobs];
   Mass[kE] = MASS_E;
   Mass[kP] = MASS_PIP;
   Mass[kM] = MASS_PIM;

   clas12::vector4 vecObs[Nobs]; // observable 4-momentum
   Float_t E[Nobs]; // energy 
   Float_t Emax[Nobs]; // maximum energy 
   Float_t P[Nobs][3]; // momentum [x,y,z]
   Int_t Idx[Nobs]; // index of observable in particleList
   Float_t Pt[Nobs]; // transverse momentum
   Float_t Phi[Nobs]; // azimuthal angle
   Float_t Theta[Nobs]; // scattering angle
   Float_t Eta[Nobs]; // pseudorapidity

   enum xyz {kX,kY,kZ};


   Bool_t foundAllObservables;


   
   while(reader.next()==true) {
     b.resume();

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
         if(particleList.getPid(i) == PID[o]) {
           particleList.getVector4(i,vecObs[o],Mass[o]);
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
     for(int o=0; o<Nobs; o++) if(Idx[o]<0) foundAllObservables=false;
     if(!foundAllObservables) continue;
     if(debug) printf("-------\n");


     // loop through observables, calculating kinematics
     for(int o=0; o<Nobs; o++) { 
       // -- set E to Emax, and set 4-momenta
       E[o] = Emax[o];
       // -- set observable momenta
       particleList.getVector4(Idx[o],vecObs[o],Mass[o]); // set 4-momentum

       P[o][kX] = particleList.getPx(Idx[o]);
       P[o][kY] = particleList.getPy(Idx[o]);
       P[o][kZ] = particleList.getPz(Idx[o]);

       //P[o][kX] = vecObs[o].vect().x(); // faster way to get momenta ??
       //P[o][kY] = vecObs[o].vect().y();
       //P[o][kZ] = vecObs[o].vect().z();
       
       Pt[o] = TMath::Hypot( P[o][kX], P[o][kY] ); // transverse momentum
       // -- compute angles
       Phi[o] = TMath::ATan2( P[o][kY], P[o][kX] );
       Theta[o] = TMath::ACos(Pt[o]);
       Eta[o] = -1 * TMath::Log( TMath::Tan( Theta[o]/2.0 ) ); 
       if(debug) printf("[%s]\tindex=%d  energy=%.3f\n",Name[o].Data(),Idx[o],E[o]);
     };


     // calculate DIS kinematics
     vecW = vecBeam + vecTarget - vecObs[kE];
     vecQ = vecBeam - vecObs[kE];
     Q2 = -1*vecQ.m2();
     W = vecW.m();
     //Nu = vecTarget.e() * vecQ.e() / MASS_P;
     Nu = vecBeam.e() - vecObs[kE].e();
     X = Q2 / ( 2 * MASS_P * Nu );

     if(debug) printf("--> W = %f\n",W);


     tree->Fill();

     b.pause();
   };
   reader.showBenchmark();
   printf(" time spend on analysis = %8.5f (ms) events = %8d\n",
   b.getTime()*1e-9,
      b.getCounter());
   tree->Write("tree");
};
//### END OF GENERATED CODE
