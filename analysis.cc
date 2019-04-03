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
static const int PID_E = 11;
static const int PID_P = 2212;
static const int PID_N = 2112;
static const int PID_PIP = 211;
static const int PID_PIM = -211;
static const int PID_PI0 = 111;
static const int PID_KP = 321;
static const int PID_KM = -321;
static const int PID_PHOTON = 22;

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


int main(int argc, char** argv) {

   // ARGUMENTS
   TString infileN,outfileN;
   if(argc>1) infileN = TString(argv[1]);
   else {
     printf("USAGE: %s [hipo file]\n",argv[0]);
     exit(0);
   };
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


   clas12::vector4 beam( 0.0, 0.0, 10.5, 10.5   );
   clas12::vector4 target( 0.0, 0.0,  0.0,  MASS_P );
   clas12::vector4 electron;
   clas12::vector4 vecW;
   clas12::vector4 vecQ;

   hipo::benchmark b;

   Float_t eleE,eleEmax;
   Int_t idxE;
   
   while(reader.next()==true) {
     b.resume();

     int size = particleList.getSize();

     // search for highest-energy electron
     eleEmax = 0;
     idxE = -1;
     for(int i=0; i<size; i++) {
       if(particleList.getPid(i) == 11) {
         particleList.getVector4(i,electron,MASS_E);
         eleE = electron.e();
         if(eleE > eleEmax) {
           eleEmax = eleE;
           idxE = i;
         };
       };
     };
     if(idxE<0) continue;


     // calculate DIS kinematics
     particleList.getVector4(idxE,electron,MASS_E);
     vecW = beam + target - electron;
     vecQ = beam - electron;
     Q2 = -1*vecQ.m2();
     W = vecW.m();
     //Nu = target.e() * vecQ.e() / MASS_P;
     Nu = beam.e() - electron.e();
     X = Q2 / ( 2 * MASS_P * Nu );


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
