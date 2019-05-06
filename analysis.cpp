#include <cstdlib>
#include <iostream>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRegexp.h"

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
   Bool_t batchMode = false;
   if(argc<=1) {
     printf("USAGE: %s [hipo file]\n",argv[0]);
     exit(0);
   };
   if(argc>1) infileN = TString(argv[1]);
   if(argc>2) batchMode = true;


   // if true, print more to stdout
   bool debug = 0;

   // set output file name
   TString outfileN;
   if(batchMode) outfileN = "outroot.root";
   else {
     outfileN = infileN;
     outfileN(TRegexp("^.*/")) = "outroot/";
     outfileN(TRegexp("hipo$")) = "root";
   };
   printf("outfileN = %s\n",outfileN.Data());


   // load libs
   DIS * disEv = new DIS();
   Trajectory * hadron[2];
   Dihadron * dih = new Dihadron(); dih->useBreit = false;
   Dihadron * dihBr = new Dihadron(); dihBr->useBreit = true;

   hadron[hP] = new Trajectory(kPip);
   hadron[hM] = new Trajectory(kPim);

   hipo::benchmark bench;


   // debugging flags
   disEv->debug = 0;
   dih->debug = 0;
   dihBr->debug = 0;


   // define tree
   TFile * outfile = new TFile(outfileN,"RECREATE");
   TTree * tree = new TTree();

   // DIS kinematics branches
   tree->Branch("W",&(disEv->W),"W/F");
   tree->Branch("Q2",&(disEv->Q2),"Q2/F");
   tree->Branch("Nu",&(disEv->Nu),"Nu/F");
   tree->Branch("x",&(disEv->x),"x/F");
   tree->Branch("y",&(disEv->y),"y/F");

   // hadron branches
   Float_t hadE[2]; // [enum plus_minus (0=+, 1=-)]
   Float_t hadP[2];
   Float_t hadPt[2];
   Float_t hadEta[2];
   Float_t hadPhi[2];
   tree->Branch("hadE",hadE,"hadE[2]/F");
   tree->Branch("hadP",hadP,"hadP[2]/F");
   tree->Branch("hadPt",hadPt,"hadPt[2]/F");
   tree->Branch("hadEta",hadEta,"hadEta[2]/F");
   tree->Branch("hadPhi",hadPhi,"hadPhi[2]/F");

   // dihadron branches
   tree->Branch("Mh",&(dih->Mh),"Mh/F");
   tree->Branch("Mmiss",&(dih->Mmiss),"Mmiss/F");
   tree->Branch("Z",dih->z,"Z[2]/F");
   tree->Branch("Zpair",&(dih->zpair),"Zpair/F");
   tree->Branch("xF",&(dih->xF),"xF/F");
   tree->Branch("alpha",&(dih->alpha),"alpha/F");
   tree->Branch("Ph",&(dih->PhMag),"Ph/F");
   tree->Branch("PhPerp",&(dih->PhPerpMag),"PhPerp/F");
   tree->Branch("R",&(dih->RMag),"R/F");
   tree->Branch("RPerp",&(dih->RPerpMag),"RPerp/F");
   tree->Branch("RT",&(dih->RTMag),"RT/F");
   tree->Branch("PhiH",&(dih->PhiH),"PhiH/F");
   tree->Branch("PhiR",&(dih->PhiR),"PhiR/F");
   // -- phiR angle test branches
   tree->Branch("PhiRq",&(dih->PhiRq),"PhiRq/F"); // via R_perp
   tree->Branch("PhiRp",&(dih->PhiRp),"PhiRp/F"); // via R_T
   tree->Branch("PhiRp_r",&(dih->PhiRp_r),"PhiRp_r/F"); // via R_T (frame-dependent)

   // breit frame dihadron branches
   tree->Branch("b_alpha",&(dihBr->alpha),"b_alpha/F");
   tree->Branch("b_Ph",&(dihBr->PhMag),"b_Ph/F");
   tree->Branch("b_PhPerp",&(dihBr->PhPerpMag),"b_PhPerp/F");
   tree->Branch("b_R",&(dihBr->RMag),"b_R/F");
   tree->Branch("b_RPerp",&(dihBr->RPerpMag),"b_RPerp/F");
   tree->Branch("b_RT",&(dihBr->RTMag),"b_RT/F");
   tree->Branch("b_PhiH",&(dihBr->PhiH),"b_PhiH/F");
   tree->Branch("b_PhiR",&(dihBr->PhiR),"b_PhiR/F");
   // -- phiR angle test branches
   tree->Branch("b_PhiRq",&(dihBr->PhiRq),"b_PhiRq/F"); // via R_perp
   tree->Branch("b_PhiRp",&(dihBr->PhiRp),"b_PhiRp/F"); // via R_T
   tree->Branch("b_PhiRp_r",&(dihBr->PhiRp_r),"b_PhiRp_r/F"); // via R_T (frame-dependent)

   // define reader and particle list
   hipo::reader reader;
   reader.open(infileN.Data());
   clas12::particle particleList("REC::Particle",reader);


   // define additional banks
   hipo::bank evBank("REC::Event",reader);
   hipo::bank configBank("RUN::config",reader);

   // Get node numbers of this bank
   //
   // NOTE: Clas12Tool Hipo/bank::getEntryOrder is protected,
   // and as far as I know, this is how to access specific bank 
   // nodes; my minimal workaround was to add the following *public* 
   // method to Hipo/bank.h:
   // 
   // int getn(const char *e) { return getEntryOrder(e); };
   // 
   Int_t o_evnum = evBank.getn("NEVENT"); // event #
   Int_t o_runnum = evBank.getn("NRUN"); // run #
   Int_t o_helicity = evBank.getn("Helic"); // e- helicity
   Int_t o_torus = configBank.getn("torus"); // torus in/outbending
   Int_t o_triggerBits = configBank.getn("trigger"); // trigger bits

   Int_t evnum,runnum;
   Int_t helicity;
   Float_t torus;
   Long64_t triggerBits;

   tree->Branch("runnum",&runnum,"runnum/I");
   tree->Branch("evnum",&evnum,"evnum/I");
   tree->Branch("helicity",&helicity,"helicity/I");
   tree->Branch("torus",&torus,"torus/F");
   tree->Branch("triggerBits",&triggerBits,"triggerBits/L");



   // define observable variables
   Float_t E[nParticles]; // energy 
   Float_t Emax[nParticles]; // maximum energy 
   Int_t Idx[nParticles]; // index of observable in particleList
   Float_t P[nParticles][3]; // momentum [x,y,z]
   enum xyz {dX,dY,dZ};
   clas12::vector4 vecObs; // observable 4-momentum
   Int_t oCur,pidCur;


   Bool_t foundAllObservables;

   Int_t evCount = 0;

   printf("begin event loop...");
   while(reader.next()==true) {
     bench.resume();


     // read event-level banks
     evnum = evBank.getInt(o_evnum,0);
     runnum = evBank.getInt(o_runnum,0);
     helicity = evBank.getInt(o_helicity,0);
     torus = configBank.getFloat(o_torus,0);
     triggerBits = configBank.getLong(o_triggerBits,0);


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
       hadE[h] = (hadron[h]->Vec).E();
       hadP[h] = (hadron[h]->Vec).P();
       hadPt[h] = (hadron[h]->Vec).Pt();
       hadEta[h] = (hadron[h]->Vec).Eta();
       hadPhi[h] = (hadron[h]->Vec).Phi();
       if(debug) {
         printf("[+] %s 4-momentum:\n",(hadron[h]->Title()).Data());
         (hadron[h]->Vec).Print();
       };
     };


     // set dihadron momenta
     dih->SetEvent(hadron[hP],hadron[hM],disEv);
     dihBr->SetEvent(hadron[hP],hadron[hM],disEv);


       


     tree->Fill();
     evCount++;
     if(evCount%100==0) 
       printf("[---] %d dihadron events found\n",evCount);

     bench.pause();
   };
   reader.showBenchmark();
   printf(" time spend on analysis = %8.5f (ms) events = %8d\n",
     bench.getTime()*1e-9,
     bench.getCounter());
   tree->Write("tree");
   printf("\nDIHADRON EVENT COUNT: %d\n\n",evCount);

   outfile->Close();
};
