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
   bool debug = 1;
   bool debugSort = 1;

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
   Dihadron * dih = new Dihadron(); dih->useBreit = false;
   Dihadron * dihBr = new Dihadron(); dihBr->useBreit = true;

   // define trajectories
   Trajectory * had[2]; // for dihadron pair; points to traj instances
   Trajectory * ele; // for electron
   Trajectory * photon[2]; // for diphoton (for pi0)

   // set up structure for sorting particles by E
   //   each particle has a trajectory pointer array,
   //   which will be sorted by energy later
   const Int_t maxTraj = 40;
   TObjArray * trajArrUS[nParticles]; // unsorted array
   TObjArray * trajArr[nParticles]; // sorted array
   Trajectory * traj[nParticles][maxTraj]; // array of trajectory pointers
   Trajectory * tr; // set to a pointer in traj array
   Int_t trajIdx[nParticles][maxTraj]; // sort index for each particle
   Float_t trajE[nParticles][maxTraj]; // sorted energy for each particle
   Int_t trajCnt[nParticles]; // number of particles for each type of observable
   
   for(int o=0; o<nParticles; o++) {
     trajArrUS[o] = new TObjArray();
     trajArr[o] = new TObjArray();
     trajCnt[o] = 0;
     for(int t=0; t<maxTraj; t++) {
       traj[o][t] = new Trajectory(o);
     };
   };
       

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

   // miscellaneous branches for classifying the type of observables
   Int_t pairType;
   TString particleCntStr = Form("particleCnt[%d]/I",nParticles);
   tree->Branch("pairType",&pairType,"pairType/I"); // dihadron pair type (see Constants.h)
   tree->Branch("particleCnt",trajCnt,particleCntStr); // number of particles


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
   tree->Branch("PhEta",&(dih->PhEta),"PhEta/F");
   tree->Branch("PhPhi",&(dih->PhPhi),"PhPhi/F");
   tree->Branch("R",&(dih->RMag),"R/F");
   tree->Branch("RPerp",&(dih->RPerpMag),"RPerp/F");
   tree->Branch("RT",&(dih->RTMag),"RT/F");
   tree->Branch("PhiH",&(dih->PhiH),"PhiH/F");
   // -- phiR angles
   tree->Branch("PhiRq",&(dih->PhiRq),"PhiRq/F"); // via R_perp
   tree->Branch("PhiRp",&(dih->PhiRp),"PhiRp/F"); // via R_T
   tree->Branch("PhiRp_r",&(dih->PhiRp_r),"PhiRp_r/F"); // via R_T (frame-dependent)
   tree->Branch("PhiRp_g",&(dih->PhiRp_g),"PhiRp_g/F"); // via eq. 9 in 1408.5721

   // breit frame dihadron branches
   tree->Branch("b_alpha",&(dihBr->alpha),"b_alpha/F");
   tree->Branch("b_Ph",&(dihBr->PhMag),"b_Ph/F");
   tree->Branch("b_PhPerp",&(dihBr->PhPerpMag),"b_PhPerp/F");
   tree->Branch("b_PhEta",&(dihBr->PhEta),"b_PhEta/F");
   tree->Branch("b_PhPhi",&(dihBr->PhPhi),"b_PhPhi/F");
   tree->Branch("b_R",&(dihBr->RMag),"b_R/F");
   tree->Branch("b_RPerp",&(dihBr->RPerpMag),"b_RPerp/F");
   tree->Branch("b_RT",&(dihBr->RTMag),"b_RT/F");
   tree->Branch("b_PhiH",&(dihBr->PhiH),"b_PhiH/F");
   // -- phiR angles
   tree->Branch("b_PhiRq",&(dihBr->PhiRq),"b_PhiRq/F"); // via R_perp
   tree->Branch("b_PhiRp",&(dihBr->PhiRp),"b_PhiRp/F"); // via R_T
   tree->Branch("b_PhiRp_r",&(dihBr->PhiRp_r),"b_PhiRp_r/F"); // via R_T (frame-dependent)
   tree->Branch("b_PhiRp_g",&(dihBr->PhiRp_g),"b_PhiRp_g/F"); // via eq. 9 in 1408.5721

   // define reader and particle list
   hipo::reader reader;
   reader.open(infileN.Data());
   clas12::particle particleList("REC::Particle",reader);
   //hipo::bank particleBank("REC::Particle",reader); // (only used for printout for debug)


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
   /*
   Float_t E[nParticles]; // energy 
   Float_t Emax[nParticles]; // maximum energy 
   Int_t Idx[nParticles]; // index of observable in particleList
   */
   clas12::vector4 vecObs; // observable 4-momentum
   Int_t oCur,pidCur;


   Bool_t foundAllObservables[nPairType];

   Int_t evCount = 0;
   Int_t pairCount[nPairType];
   for(int pp=0; pp<nPairType; pp++) pairCount[pp]=0;

   printf("begin event loop...");
   while(reader.next()==true) {
     bench.resume();


     // read event-level banks
     evnum = evBank.getInt(o_evnum,0);
     runnum = evBank.getInt(o_runnum,0);
     helicity = evBank.getInt(o_helicity,0);
     torus = configBank.getFloat(o_torus,0);
     triggerBits = configBank.getLong(o_triggerBits,0);

     //if(debugSort) particleBank.show();
     if(debugSort) { for(int l=0; l<30; l++) printf("-"); printf("\n"); };


     // search for observables, sorting each kind by energy
     // ---------------------------------------------------
     //
     // -- reset trajectory-sorting data structures
     for(int o=0; o<nParticles; o++) {
       trajCnt[o] = 0;
       trajArrUS[o]->Clear();
       trajArr[o]->Clear();
       for(int t=0; t<maxTraj; t++) {
         trajIdx[o][t] = -1;
         trajE[o][t] = -1;
       };
     };
     
     // -- if it's an observable that we want, add its trajectory
     //    to the unsorted trajectory array and its energy to list
     //    of energies for this observable
     for(int i=0; i<particleList.getSize(); i++) {

       pidCur = particleList.getPid(i);
       if      (pidCur == PartPID(kE)   )  oCur=kE;
       else if (pidCur == PartPID(kPip) )  oCur=kPip;
       else if (pidCur == PartPID(kPim) )  oCur=kPim;
       else if (pidCur == PartPID(kPi0) )  oCur=kPi0;
       else if (pidCur == PartPID(kPhoton) ) oCur=kPhoton;
       else oCur=-10000;

       if(debugSort) printf(" pid=%d  oCur=%d\n",pidCur,oCur);

       if(oCur>-10000) {

         particleList.getVector4(i,vecObs,PartMass(oCur));

         // set Trajectory pointer tr to proper allocated Trajectory instance;
         // if there are more instances of this observable than we allocated for,
         // we could allocate more memory for the extra ones (which could cause
         // a memory leak or eventual segfault if maxTraj is not high enough...),
         // or ignore the extra particles
         if(trajCnt[oCur]<maxTraj) {
           tr = traj[oCur][trajCnt[oCur]]; // set tr to allocated Trajectory instance
         } else {
           fprintf(stderr,"WARNING: more than maxTraj observables of type %s found\n",
             PartName(oCur).Data());
           //tr = new Trajectory(oCur); // allocate more memory
           continue; // ignore this particle 
         };

         // add this Trajectory to the unsorted array
         tr->SetMomentum(
           particleList.getPx(i),
           particleList.getPy(i),
           particleList.getPz(i)
         );
         trajArrUS[oCur]->AddLast(tr);

         // add this trajectory's energy to the energy array
         trajE[oCur][trajCnt[oCur]] = vecObs.e();

         // increment the trajectory counter
         trajCnt[oCur]++;

         if(debugSort) {
           printf("found %s (pid=%d):\n",PartName(oCur).Data(),pidCur);
           tr->Vec.Print();
         };
       };
     };
     if(debugSort) printf(">>\n");

     // -- sort each list of each type of observable by E; trajArr will 
     //    be the sorted trajectory array and the 0th entry will be the 
     //    highest energy trajectory 
     for(int o=0; o<nParticles; o++) { 
       if(trajCnt[o]>0) {

         // sort the energy array
         TMath::Sort(trajCnt[o],trajE[o],trajIdx[o]);

         // sort the trajectory array according to the sorted energy array
         for(int t=0; t<trajCnt[o]; t++) {
           tr = (Trajectory*) trajArrUS[o]->At(trajIdx[o][t]);
           trajArr[o]->AddLast(tr);
         };

         if(debugSort) {
           printf("sorted list: %s cnt=%d\n",PartName(o).Data(),trajCnt[o]);
           for(int t=0; t<trajCnt[o]; t++)
             ((Trajectory*)(trajArr[o]->At(t)))->Vec.Print();
         };
       };

     };
     //
     // ---------------------------------------------------



     // check if we found all the observables for a given pairType
     foundAllObservables[pairPM] = trajCnt[kE]>0 && trajCnt[kPip]>0 && trajCnt[kPim]>0;
     foundAllObservables[pairP0] = trajCnt[kE]>0 && trajCnt[kPip]>0 && trajCnt[kPi0]>0;
     foundAllObservables[pair0M] = trajCnt[kE]>0 && trajCnt[kPi0]>0 && trajCnt[kPim]>0;

     
     // loop through pair types
     for(int pp=0; pp<nPairType; pp++) {

       if(foundAllObservables[pp]) {

         pairType = pp;

         if(debug) printf(">>> BEGIN DIHADRON EVENT %s%s\n",
           PMsym(pairType,hP).Data(),PMsym(pairType,hM).Data());

         // compute DIS kinematics
         ele = (Trajectory*) trajArr[kE]->At(0); // select highest-E electron
         disEv->SetElectron(ele);
         disEv->Analyse();


         // set hadron kinematics
         for(int h=0; h<2; h++) {

           // select highest-E hadron
           had[h] = (Trajectory*) trajArr[PMidx(pairType,h)]->At(0); 

           hadE[h] = (had[h]->Vec).E();
           hadP[h] = (had[h]->Vec).P();
           hadPt[h] = (had[h]->Vec).Pt();
           hadEta[h] = (had[h]->Vec).Eta();
           hadPhi[h] = (had[h]->Vec).Phi();
           if(debug) {
             printf("[+] %s 4-momentum:\n",(had[h]->Title()).Data());
             (had[h]->Vec).Print();
           };
         };


         // compute dihadron kinematics
         dih->SetEvent(had[hP],had[hM],disEv);
         dihBr->SetEvent(had[hP],had[hM],disEv);


         // fill tree
         tree->Fill();

         // increment event counter
         evCount++;
         pairCount[pairType]++;
         if(evCount%100==0) printf("[---] %d dihadron events found\n",evCount);

       }; // eo if foundAllObservables for this pairType
     }; // eo pairType loop

     bench.pause();
   };
   reader.showBenchmark();
   printf(" time spend on analysis = %8.5f (ms) events = %8d\n",
     bench.getTime()*1e-9,
     bench.getCounter());
   tree->Write("tree");
   printf("\nDIHADRON EVENT COUNT: %d\n",evCount);
   for(int pp=0; pp<nPairType; pp++) 
     printf(" %s : %d\n",pairName(pp).Data(),pairCount[pp]);
   printf("\n");


   outfile->Close();
};
