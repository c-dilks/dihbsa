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
#include "TLorentzVector.h"


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


int main(int argc, char** argv) {

#ifdef HIPO_VERSION
   printf("%s compiled with HIPO_VERSION = %d\n",argv[0],HIPO_VERSION);
#else
   fprintf(stderr,"ERROR: HIPO_VERSION preprocessor macro undefined\n");
   exit(0);
#endif

   
   // ARGUMENTS
   TString infileN;
   Bool_t batchMode = false;
   if(argc<=1) {
     printf("USAGE: %s [hipo file]\n",argv[0]);
     exit(0);
   };
   if(argc>1) infileN = TString(argv[1]);
   if(argc>2) batchMode = true;


   // debugging flags
   bool debug = 0; // general debugging statements
   bool debugPHPsort = 0; // sorting photon pairs by energy





   // set output file
   TString outfileN;
   if(batchMode) outfileN = "outroot.root";
   else {
     outfileN = infileN;
     outfileN(TRegexp("^.*/")) = "outroot/";
     outfileN(TRegexp("hipo$")) = "root";
   };
   printf("outfileN = %s\n",outfileN.Data());
   TFile outfile(outfileN,"RECREATE");


   // load libs
   DIS * disEv = new DIS();
   Dihadron * dih = new Dihadron(); dih->useBreit = false;
   Dihadron * dihBr = new Dihadron(); dihBr->useBreit = true;
   Diphoton * diPhot[2];
   for(int h=0; h<2; h++) diPhot[h] = new Diphoton();

   // define trajectories
   Trajectory * had[2]; // for dihadron pair; points to traj instances
   Trajectory * ele; // for electron
   Trajectory * phot[2]; // for diphoton (for pi0)

   // set up structure for sorting particles by E
   //   each particle has a trajectory pointer array,
   //   which will be sorted by energy later
   const Int_t maxTraj = 20;
   TObjArray * trajArrUS[nParticles]; // unsorted array
   TObjArray * trajArr[nParticles]; // sorted array
   Trajectory * traj[nParticles][maxTraj]; // array of trajectory pointers
   Trajectory * tr; // set to a pointer in traj array
   Int_t trajSortIdx[nParticles][maxTraj]; // sort index for each particle
   Float_t trajE[nParticles][maxTraj]; // sorted energy for each particle
   Int_t trajCnt[nParticles]; // number of particles for each type of observable
   Int_t nanCnt[nParticles]; // counts number of particles which have NaN 4-momenta
   Int_t gtMaxTrajCnt[nParticles]; // counts number of particles that have multiplicity
                                   // greater than maxTraj
   
   for(int p=0; p<nParticles; p++) {
     trajArrUS[p] = new TObjArray();
     trajArr[p] = new TObjArray();
     trajCnt[p] = 0;
     nanCnt[p] = 0;
     gtMaxTrajCnt[p] = 0;
     for(int t=0; t<maxTraj; t++) {
       traj[p][t] = new Trajectory(p);
     };
   };
       

   hipo::benchmark bench;


   // debugging flags
   disEv->debug = 0;
   dih->debug = 0;
   dihBr->debug = 0;


   // define tree
   TTree * tree = new TTree("tree","tree");

   // DIS kinematics branches
   tree->Branch("W",&(disEv->W),"W/F");
   tree->Branch("Q2",&(disEv->Q2),"Q2/F");
   tree->Branch("Nu",&(disEv->Nu),"Nu/F");
   tree->Branch("x",&(disEv->x),"x/F");
   tree->Branch("y",&(disEv->y),"y/F");
   
   // electron kinematics branches
   tree->Branch("eleE",&(disEv->eleE),"eleE/F");
   tree->Branch("eleP",&(disEv->eleP),"eleP/F");
   tree->Branch("elePt",&(disEv->elePt),"elePt/F");
   tree->Branch("eleEta",&(disEv->eleEta),"eleEta/F");
   tree->Branch("elePhi",&(disEv->elePhi),"elePhi/F");


   // miscellaneous branches for classifying the type of observables
   Int_t particleCntAll;
   TString particleCntStr = Form("particleCnt[%d]/I",nParticles);
   tree->Branch("particleCnt",trajCnt,particleCntStr); // number of each type of particle
                                                     // that we consider, as defined in
                                                     // Constants.h
   tree->Branch("particleCntAll",&particleCntAll,"particleCntAll/I"); // total number
                                                     // of particles in the event (note:
                                                     // includes particles not in
                                                     // Constants.h)


   // hadron branches // [ordered such that first hadron has higher charge than second]
   Int_t pairType;
   Int_t hadIdx[2]; // particle Idx of each hadron in the pair
   Float_t hadE[2];
   Float_t hadP[2];
   Float_t hadPt[2];
   Float_t hadEta[2];
   Float_t hadPhi[2];
   tree->Branch("pairType",&pairType,"pairType/I");
   tree->Branch("hadIdx",hadIdx,"hadIdx[2]/I");
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
   tree->Branch("theta",&(dih->theta),"theta/F");

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
   /*
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
   */

   // pi0 (diphoton) branches
   // -- they are arrays of length diphCnt; there are three cases:
   //    - diphCnt=0: no diphotons in the dihadron
   //    - diphCnt=1: one dihadron hadrons is a diphoton (pi0 or BG)
   //    - diphCnt=2: both dihadron hadrons are diphotons (pi0-pi0, pi0-BG, or BG-BG)
   Int_t diphCnt,diphCnt_tr;
   Float_t diphPhotE[2][2], diphPhotPt[2][2], diphPhotEta[2][2], diphPhotPhi[2][2];
   Float_t diphE[2], diphZ[2], diphPt[2], diphM[2], diphAlpha[2], diphEta[2], diphPhi[2];
   tree->Branch("diphCnt",&diphCnt_tr,"diphCnt/I"); // number of diphotons {0,1,2}
   tree->Branch("diphPhotE",   diphPhotE,   "diphPhotE[diphCnt][2]/F"); // photon energy
   tree->Branch("diphPhotPt",  diphPhotPt,  "diphPhotPt[diphCnt][2]/F"); // photon pT
   tree->Branch("diphPhotEta", diphPhotEta, "diphPhotEta[diphCnt][2]/F"); // photon eta
   tree->Branch("diphPhotPhi", diphPhotPhi, "diphPhotPhi[diphCnt][2]/F"); // photon phi
   tree->Branch("diphE",     diphE,     "diphE[diphCnt]/F"); // energy
   tree->Branch("diphZ",     diphZ,     "diphZ[diphCnt]/F"); // energy imbalance
   tree->Branch("diphPt",    diphPt,    "diphPt[diphCnt]/F"); // pT
   tree->Branch("diphM",     diphM,     "diphM[diphCnt]/F"); // invariant mass
   tree->Branch("diphAlpha", diphAlpha, "diphAlpha[diphCnt]/F"); // opening angle
   tree->Branch("diphEta",   diphEta,   "diphEta[diphCnt]/F"); // eta
   tree->Branch("diphPhi",   diphPhi,   "diphPhi[diphCnt]/F"); // pT

   // miscellaneous event-header branches
   Int_t evnum,runnum;
   Int_t helicity;
   Float_t torus,solenoid;
   Long64_t triggerBits;
   tree->Branch("runnum",&runnum,"runnum/I");
   tree->Branch("evnum",&evnum,"evnum/I");
   tree->Branch("helicity",&helicity,"helicity/I");
   tree->Branch("torus",&torus,"torus/F");
   tree->Branch("solenoid",&solenoid,"solenoid/F");
   tree->Branch("triggerBits",&triggerBits,"triggerBits/L");




   // define HIPO file reader and banks
#if HIPO_VERSION == 3
   // reader
   hipo::reader reader; // HIPO3
   reader.open(infileN.Data());

   // particle bank
   // -- for HIPO3, need to use general bank rather than clas12::particle
   //    (see log 14.6.19 for details; seems HIPO3's clas12::particle is a bit broken)
   hipo::bank particleBank("REC::Particle",reader);  // HIPO3
   Int_t o_pid = particleBank.getn("pid");
   Int_t o_px = particleBank.getn("px");
   Int_t o_py = particleBank.getn("py");
   Int_t o_pz = particleBank.getn("pz");
   /*
   clas12::particle particleBank("REC::Particle",reader); // see log 14.6.19
   particleBank.init("REC::Particle",reader); // init needs to be called in HIPO3 version
   */

   // event and run config banks
   // -- In Hipo3 Clas12Tool code, to access contents of bank entries, the entry order
   //    number within the bank is needed, which can only be obtained by
   //    Clas12Tool/Hipo3/bank::getEntryOrder, a protected method.
   //
   // -- simple workaround: add the following public method to Hipo3/bank.h:
   //
   //   int getn(const char *e) { return getEntryOrder(e); };
   //
   // -- in Hipo4/bank.h, bank entries are now accessible by name
   //
   hipo::bank configBank("RUN::config",reader); // HIPO3
   hipo::bank evBank("REC::Event",reader); // HIPO3
   Int_t o_torus = configBank.getn("torus"); // torus in/outbending
   Int_t o_triggerBits = configBank.getn("trigger"); // trigger bits
   Int_t o_evnum = evBank.getn("NEVENT"); // event #
   Int_t o_runnum = evBank.getn("NRUN"); // run #
   Int_t o_helicity = evBank.getn("Helic"); // e- helicity

#elif HIPO_VERSION == 4
   clas12::clas12reader reader(infileN.Data()); // HIPO4
   printf("BEGIN TEST DICTIONARY READ\n");
   hipo::dictionary factory;
   reader.getReader().readDictionary(factory);
   //factory.show();

   hipo::event readerEvent;
   hipo::bank mcLund(factory.getSchema("MC::Lund"));
#endif

   


   // define observable variables
   TLorentzVector vecObs; // observable 4-momentum
   Float_t vecObsPx,vecObsPy,vecObsPz;
   Int_t pidCur,pIdx;
   Int_t i1,i2;

   Bool_t foundObservablePair;


   // photon pair ("php") variables
   Int_t phpCnt; // number of photon pairs
   const Int_t phpCntMax = 15; // max number of photon pairs considered
   Int_t phpIdx[phpCntMax][2]; // index of each photon in the pair
   Float_t phpE[phpCntMax]; // energy of photon pair (to be sorted)
   Int_t phpSortIdx[phpCntMax]; // sorted pair index
   Int_t phpSI[2]; // index of photon from the sorted pair list
   Bool_t photonUsed[maxTraj];

   Diphoton * diPhotTmp[phpCntMax]; // used to compute kinematics for photon pairs
   for(int php=0; php<phpCntMax; php++) diPhotTmp[php] = new Diphoton();


   // event and pair counters
   Int_t evCount = 0;



   // ----------------------------------------------------
   // EVENT LOOP
   // ----------------------------------------------------
   printf("begin event loop...\n");

   while(reader.next()==true) {
     //if(evCount>1e6) { fprintf(stderr,"BREAKING LOOP HERE!!!\n"); break; };
     bench.resume();


     // reset branches
     disEv->ResetVars();
     dih->ResetVars();
     dihBr->ResetVars();
     for(int h=0; h<2; h++) {
       hadIdx[h] = -10000;
       hadE[h] = -10000;
       hadP[h] = -10000;
       hadPt[h] = -10000;
       hadEta[h] = -10000;
       hadPhi[h] = -10000;
     };
     pairType = -10000;
     for(int h=0; h<2; h++) diPhot[h]->ResetVars();
     for(int php=0; php<phpCntMax; php++) diPhotTmp[php]->ResetVars();
     diphCnt = 0;



     // read event-header stuff
     //    TODO: I'm guessing the same HIPO4 accessors now exist for the HIPO3 version
     //    of Clas12Tool, so we may be able to get rid of this preprocessor conditional
#if HIPO_VERSION == 3
     evnum = evBank.getInt(o_evnum,0); // -->tree
     runnum = evBank.getInt(o_runnum,0); // -->tree
     helicity = evBank.getInt(o_helicity,0); // -->tree
     triggerBits = configBank.getLong(o_triggerBits,0); // -->tree
     torus = configBank.getFloat(o_torus,0); // -->tree
     torus = -10000;
     solenoid = -10000;
#elif HIPO_VERSION == 4
     evnum = reader.runconfig()->getEvent(); // -->tree
     runnum = reader.runconfig()->getRun(); // -->tree
     helicity = reader.event()->getHelicity(); // -->tree
     triggerBits = reader.runconfig()->getTrigger(); // -->tree
     torus = reader.runconfig()->getTorus(); // -->tree
     solenoid = reader.runconfig()->getSolenoid(); // -->tree
     //printf("schema name = %s\n",reader.head()->getSchema().getName().data());
#endif


     if(debug) Tools::PrintSeparator(30);
     //if(debug) particleBank.show(); // causes segfault!


     // reset trajectory-sorting data structures "trajArr" and "trajE"
     for(int p=0; p<nParticles; p++) {
       trajCnt[p] = 0;
       trajArrUS[p]->Clear();
       trajArr[p]->Clear();
       for(int t=0; t<maxTraj; t++) {
         trajSortIdx[p][t] = -1;
         trajE[p][t] = -1;
       };
     };

     
     // ---------------------------------------------------
     // PARTICLE LOOP
     // ---------------------------------------------------
     // -- read in each particle and put them into trajArr, which will be sorted
     //    afterward
     // -- the way we loop through particles differs between HIPO versions (for now...)
     //    TODO: I'm guessing the same HIPO4 particle momentum accessors now exist for
     //    the HIPO3 version of Clas12Tool, so we may be able to get rid of this
     //    preprocessor conditional
     
#if HIPO_VERSION == 3
     particleCntAll = particleBank.getSize(); // -->tree
     if(debug) printf("particleBank.getSize() = %d\n",particleBank.getSize());
     for(int i=0; i<particleCntAll; i++) {

       // get particle PID and momentum components
       pidCur = particleBank.getInt(o_pid,i);
       vecObsPx = particleBank.getFloat(o_px,i);
       vecObsPy = particleBank.getFloat(o_py,i);
       vecObsPz = particleBank.getFloat(o_pz,i);

#elif HIPO_VERSION == 4
     /*
     // reconstructed particles
     particleCntAll = reader.getNParticles(); // -->tree
     if(debug) printf("reader.getNParticles() = %d\n",particleCntAll);
     for(auto & part : reader.getDetParticles()) {
       // get particle PID and momentum components
       pidCur = part->getPid();
       vecObsPx = part->par()->getPx();
       vecObsPy = part->par()->getPy();
       vecObsPz = part->par()->getPz();
       */
     // MC::Lund particles
     reader.getReader().read(readerEvent);
 		 readerEvent.getStructure(mcLund);
     particleCntAll = mcLund.getRows(); // -->tree
     printf("mcLund.getRows() = %d\n",particleCntAll);
     for(int rr=0; rr<particleCntAll; rr++) {
       pidCur = mcLund.getInt("pid",rr);
       vecObsPx = mcLund.getFloat("px",rr);
       vecObsPy = mcLund.getFloat("py",rr);
       vecObsPz = mcLund.getFloat("pz",rr);
#endif


       // convert PID to local particle index; if it's not defined in Constants.h, pIdx
       // will be -10000 and this particle will be ignored
       pIdx = PIDtoIdx(pidCur);
       if(pIdx==kP || pIdx==kN) pIdx=-10000; // also skip protons and neutrons
       if(debug) printf(" pid=%d  pIdx=%d\n",pidCur,pIdx);


       if(pIdx>-10000) {

         // set Trajectory pointer tr to proper allocated Trajectory instance; if there
         // are more instances of this observable than we allocated for, (this max
         // number is "maxTraj"), just ignore the additional particles and print a
         // warning to stderr
         if(trajCnt[pIdx]<maxTraj) {
           
           // make sure particle has sensible momentum (!=NaN)
           // (sometimes photons in dnp2018 skim files don't...)
           if( !isnan(vecObsPx) && !isnan(vecObsPy) && !isnan(vecObsPz) ) {

             // set tr to allocated Trajectory instance and add it to the unsorted array
             tr = traj[pIdx][trajCnt[pIdx]]; 
             vecObs.SetXYZM(vecObsPx,vecObsPy,vecObsPz,PartMass(pIdx));
             tr->SetVec(vecObs);
             trajArrUS[pIdx]->AddLast(tr);

             // add this trajectory's energy to the energy array
             trajE[pIdx][trajCnt[pIdx]] = vecObs.E();

             // increment the trajectory counter
             trajCnt[pIdx]++; // -->tree

             if(debug) {
               printf("found %s (pid=%d):\n",PartName(pIdx).Data(),pidCur);
               tr->Vec.Print();
             };
           } else {
             nanCnt[pIdx]++; // count instances of NaN particles
           };
         } else {
           //fprintf(stderr,"WARNING: more than maxTraj observables of type %s found\n",
             //PartName(pIdx).Data());
           //tr = new Trajectory(pIdx); // allocate more memory
           gtMaxTrajCnt[pIdx]++;
         };
       };
     }; // eo particle loop




     // -- sort each list of each type of observable by E; trajArr will 
     //    be the sorted trajectory array and the 0th entry will be the 
     //    highest energy trajectory 
     if(debug) printf(">>\n");
     for(int p=0; p<nParticles; p++) { 
       if(trajCnt[p]>0) {

         // sort the energy array
         TMath::Sort(trajCnt[p],trajE[p],trajSortIdx[p]);

         // fill the trajectory array according to the sorted energy array
         for(int t=0; t<trajCnt[p]; t++) {
           tr = (Trajectory*) trajArrUS[p]->At(trajSortIdx[p][t]);
           trajArr[p]->AddLast(tr);
         };

         if(debug) {
           printf("sorted list: %s cnt=%d\n",PartName(p).Data(),trajCnt[p]);
           for(int t=0; t<trajCnt[p]; t++)
             ((Trajectory*)(trajArr[p]->At(t)))->Vec.Print();
         };
       };

     };

     
     // ---------------------------------------------------
     // PHOTON PAIRING, for pi0s
     // ---------------------------------------------------
     // -- this is done by sorting all photon pairs by energy, and picking the
     //    highest-E pair that satisfies basic requirements set in Diphoton class
     //    booleans (e.g., an opening angle cut)

     // if a pi0 has been found in the particle bank, its trajectory will already be in
     // trajArr, and the highest-E pi0 found here will be the one sent into the tree;
     // the user will be warned of this happening
     if(trajCnt[kPi0] > 0) {
       fprintf(stderr,
         "WARNING WARNING WARNING WARNING WARNING: found pi0 in paricle bank!!!\n");
       // note that EventTree is not yet programmed to accept these pions...
     };


     // reset data structures for photon pair sorting
     phpCnt = 0;
     for(int php=0; php<phpCntMax; php++) {
       phpIdx[php][0] = -1;
       phpIdx[php][1] = -1;
       phpE[php] = -1;
       phpSortIdx[php] = -1;
     };
     for(int pp=0; pp<maxTraj; pp++) photonUsed[pp] = false;


     // only do pairing if there are at least 2 photons
     if( trajCnt[kPhoton] >= 2 ) {
       
       if(debugPHPsort) printf(">>>>>>\n");

       // loop through pairs of photons; allows for a maximum of phpCntMax pairs
       for(int p0=0; p0<trajCnt[kPhoton]; p0++) {
         for(int p1=p0+1; p1<trajCnt[kPhoton]; p1++) {
           if(phpCnt < phpCntMax) {

             // fill php index and energy arrays
             phpIdx[phpCnt][0] = p0;
             phpIdx[phpCnt][1] = p1;
             phot[0] = (Trajectory*) trajArr[kPhoton]->At(p0);
             phot[1] = (Trajectory*) trajArr[kPhoton]->At(p1);
             phpE[phpCnt] = (phot[0]->Vec).E() + (phot[1]->Vec).E();

             if(debugPHPsort) {
               printf("phpIdx[%d] = (%d,%d)  E = %.3f + %.3f = %.3f\n",
                 phpCnt, p0, p1, (phot[0]->Vec).E(), (phot[1]->Vec).E(), phpE[phpCnt]);
             };

             phpCnt++;

           };
         };
       };


       // sort photon pairs by E
       TMath::Sort(phpCnt,phpE,phpSortIdx);

       if(debugPHPsort) {
         printf("-----  sorted:  -----\n");
         for(int php=0; php<phpCnt; php++) {
           for(int h=0; h<2; h++) phpSI[h] = phpIdx[phpSortIdx[php]][h];
           printf("(%d,%d)  E=%.4f\n",phpSI[0],phpSI[1],phpE[phpSortIdx[php]]);
         };
       };


       // loop through sorted pair list, adding satisfactory diphotons to trajArr
       if(phpCnt>0) { 
         for(int php=0; php<phpCnt; php++) {

           // set phot Trajectories to this pair of photons
           for(int h=0; h<2; h++) {
             phpSI[h] = phpIdx[phpSortIdx[php]][h];
             phot[h] = (Trajectory*) trajArr[kPhoton]->At(phpSI[h]);
           };

           // check these photons have not already been added to diphoton pair array
           if(!photonUsed[phpSI[0]] && !photonUsed[phpSI[1]]) {

             // if the pair satisfies basic requirements, add it to trajArr; 
             diPhotTmp[php]->SetEvent(phot[0],phot[1]); // compute kinematics
             if(diPhotTmp[php]->validDiphoton) {
               for(int h=0; h<2; h++) photonUsed[phpSI[h]] = true; // mark photons used
               trajArr[kDiph]->AddLast(diPhotTmp[php]->Traj);

               // if it's the highest-E one in a dihadron with 1 or 2 diphotons, or the
               // second highest-E one in a dihadron with 2 diphotons, set diPhot
               // to this pair and increment diphCnt
               if(trajCnt[kDiph]==0 || trajCnt[kDiph]==1) {
                 diPhot[trajCnt[kDiph]]->SetEvent(phot[0],phot[1]);
                 diphCnt++;
               };

               trajCnt[kDiph]++;
             };
           };

         };
       };
     };


     // fill diphoton tree branches if we have diphotons
     if(diphCnt>2) {
       fprintf(stderr,"ERROR: diphCnt somehow greater than 2; setting to 2...\n");
       diphCnt = 2;
     };
     if(diphCnt>0) {
       for(int dp=0; dp<diphCnt; dp++) {
         for(int h=0; h<2; h++) {
           diphPhotE[dp][h] = diPhot[dp]->photE[h]; // -->tree
           diphPhotPt[dp][h] = diPhot[dp]->photPt[h]; // -->tree
           diphPhotEta[dp][h] = diPhot[dp]->photEta[h]; // -->tree
           diphPhotPhi[dp][h] = diPhot[dp]->photPhi[h]; // -->tree
         };
         diphE[dp] = diPhot[dp]->E; // -->tree
         diphZ[dp] = diPhot[dp]->Z; // -->tree
         diphPt[dp] = diPhot[dp]->Pt; // -->tree
         diphM[dp] = diPhot[dp]->M; // -->tree
         diphAlpha[dp] = diPhot[dp]->Alpha; // -->tree
         diphEta[dp] = diPhot[dp]->Eta; // -->tree
         diphPhi[dp] = diPhot[dp]->Phi; // -->tree
       };
     };


     

     // ---------------------------------------------------
     // HADRON PAIRING, and fill the tree
     // ---------------------------------------------------
     
     // first make sure there's an electron
     if(trajCnt[kE] > 0) {

       // compute DIS kinematics
       ele = (Trajectory*) trajArr[kE]->At(0); // select highest-E electron
       disEv->SetElectron(ele);
       disEv->Analyse(); // -->tree


       // look for "observable pairs" -- these are pairs that are used to form
       // form dihadrons; only the desired observables are paired (see observable_enum
       // in Constants.h)
       for(Int_t o1=0; o1<nObservables; o1++) {
         for(Int_t o2=o1; o2<nObservables; o2++) {

           foundObservablePair = false;

           // convert observable index to particle index
           i1 = OI(o1);
           i2 = OI(o2);

           if(i1==i2) {
             // if the indices are the same, and we have at least 2 of this observable,
             // set 1st hadron to the highest energy one, and 2nd to the next highest
             if( trajCnt[i1] >=2 ) {
               for(int h=0; h<2; h++) {
                 hadIdx[h] = i1; // -->tree
                 had[h] = (Trajectory*) trajArr[i1]->At(h);
               };
               foundObservablePair = true;
             };

           } else {
             // if the indices are different, set hadrons according to dihHadIdx
             // (see Constants.h; higher charge comes first, and if charges are
             // equal, higher mass comes first)
             if( trajCnt[i1] >= 1 && trajCnt[i2] >= 1 ) {
               for(int h=0; h<2; h++) {
                 hadIdx[h] = dihHadIdx(i1,i2,h); // -->tree
                 had[h] = (Trajectory*) trajArr[hadIdx[h]]->At(0);
               };
               foundObservablePair = true;
             };
           };


           // if a hadron pair was found, proceed with kinematics calculations
           // and fill the tree
           if(foundObservablePair) {

             if(debug) printf(">>> BEGIN DIHADRON EVENT %s\n",PairName(i1,i2).Data());

             // set hadron kinematics
             pairType = EncodePairType(hadIdx[qA],hadIdx[qB]); // -->tree
             for(int h=0; h<2; h++) {
               hadE[h] = (had[h]->Vec).E(); // -->tree
               hadP[h] = (had[h]->Vec).P(); // -->tree
               hadPt[h] = (had[h]->Vec).Pt(); // -->tree

               if(hadE[h]>0 && hadPt[h]>0) {
                 hadEta[h] = (had[h]->Vec).Eta(); // -->tree
                 hadPhi[h] = (had[h]->Vec).Phi(); // -->tree
               } else {
                 hadEta[h] = -10000;
                 hadPhi[h] = -10000;
               };

               if(debug) {
                 printf("[+] %s 4-momentum:\n",(had[h]->Title()).Data());
                 (had[h]->Vec).Print();
               };
             };


             // compute dihadron kinematics
             dih->SetEvent(had[qA],had[qB],disEv); // -->tree
             dihBr->SetEvent(had[qA],had[qB],disEv); // -->tree

             // set diphCnt_tr to corrected diphCnt, for filling diphoton branches
             // -- this solves a subtle case: assume we have a pi+ and two diphotons;
             // this will resolve to two dihadrons: A=(pi+,diphoton) and
             // B=(diphoton,diphoton), but diphCnt=2 for this event. We thus need to
             // force diphCnt to 1 for dihadron A, and to 2 for dihadron B
             diphCnt_tr = 0; // -->tree
             if(diphCnt>0) {
               if(hadIdx[qA]==kDiph || hadIdx[qB]==kDiph) diphCnt_tr = 1; // -->tree
               if(hadIdx[qA]==kDiph && hadIdx[qB]==kDiph) diphCnt_tr = 2; // -->tree
             };



             // fill tree
             tree->Fill();

             // increment event counter
             evCount++;
             if(evCount%100000==0) printf("[---] %d dihadron events found\n",evCount);


           }; // eo if foundObservablePair


         }; // eo for o2
       }; // eo for o1
     }; // eo if #electrons > 0


     bench.pause();
   }; // eo event loop
   //---------------------------------------- END EVENT LOOP
   

   // print benchmark
   printf(" time spend on analysis = %8.5f (s) events = %8d\n",
     bench.getTime()*1e-9,
     bench.getCounter());

   // write output tree
   printf("writing tree...\n");
   outfile.cd();
   tree->Write("tree");
   printf("tree written\n");

   // print number of NaN particles (if nonzero)
   for(int p=0; p<nParticles; p++) {
     if(nanCnt[p]>0) {
       fprintf(stderr,"WARNING: %d NaN %s 4-momenta were found (and ignored)\n",
         nanCnt[p],PartName(p).Data());
     };
   };

   // print number of particles greater than maxTrajCnt
   for(int p=0; p<nParticles; p++) {
     if(gtMaxTrajCnt[p]>0) {
       fprintf(stderr,"WARNING: %d instances of %s greater than maxTraj (ignored)\n",
         gtMaxTrajCnt[p],PartName(p).Data());
     };
   };

   // print dihadron event count
   printf("\nDIHADRON EVENT COUNT: %d\n\n",evCount);


   // close the output ROOT file
   outfile.Close();
};
