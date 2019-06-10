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
#include "Diphoton.h"



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
   bool debug = 0; // general debugging statements
   bool debugPHPsort = 0; // sorting photon pairs by energy

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
   Diphoton * diPhot = new Diphoton();

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
   TFile * outfile = new TFile(outfileN,"RECREATE");
   TTree * tree = new TTree();

   // DIS kinematics branches
   tree->Branch("W",&(disEv->W),"W/F");
   tree->Branch("Q2",&(disEv->Q2),"Q2/F");
   tree->Branch("Nu",&(disEv->Nu),"Nu/F");
   tree->Branch("x",&(disEv->x),"x/F");
   tree->Branch("y",&(disEv->y),"y/F");
   
   // electron kinematics branches
   tree->Branch("eleE",&(disEv->eleE),"eleE/F");
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

   // pi0 (diphoton) branches
   tree->Branch("photE",diPhot->photE,"photE[2]/F"); // photon energy
   tree->Branch("photPt",diPhot->photPt,"photPt[2]/F"); // photon transverse momentum
   tree->Branch("photEta",diPhot->photEta,"photEta[2]/F"); // photon pseudorapidity
   tree->Branch("photPhi",diPhot->photPhi,"photPhi[2]/F"); // photon azimuth
   tree->Branch("diphE",&(diPhot->E),"diphE/F"); // diphoton energy
   tree->Branch("diphZ",&(diPhot->Z),"diphZ/F"); // diphoton energy sharing
   tree->Branch("diphPt",&(diPhot->Pt),"diphPt/F"); // diphoton transverse momentum
   tree->Branch("diphM",&(diPhot->Mgg),"diphM/F"); // diphoton invariant mass
   tree->Branch("diphAlpha",&(diPhot->Alpha),"diphAlpha/F"); // diphoton opening angle
   tree->Branch("diphEta",&(diPhot->Eta),"diphEta/F"); // diphoton pseudorapidity
   tree->Branch("diphPhi",&(diPhot->Phi),"diphPhi/F"); // diphoton pseudorapidity


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
   clas12::vector4 vecObs; // observable 4-momentum
   Int_t pidCur,pIdx;


   Bool_t foundObservablePair;
   Int_t i1,i2;


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



   // EVENT LOOP ----------------------------------------------
   printf("begin event loop...\n");
   while(reader.next()==true) {
     bench.resume();

     // reset branches
     disEv->ResetVars();
     dih->ResetVars();
     dihBr->ResetVars();
     diPhot->ResetVars();
     for(int php=0; php<phpCntMax; php++) diPhotTmp[php]->ResetVars();
     for(int h=0; h<0; h++) {
       hadIdx[h] = -10000;
       hadE[h] = -10000;
       hadP[h] = -10000;
       hadPt[h] = -10000;
       hadEta[h] = -10000;
       hadPhi[h] = -10000;
     };
     pairType = -10000;



     // read event-level banks
     evnum = evBank.getInt(o_evnum,0); // -->tree
     runnum = evBank.getInt(o_runnum,0); // -->tree
     helicity = evBank.getInt(o_helicity,0); // -->tree
     torus = configBank.getFloat(o_torus,0); // -->tree
     triggerBits = configBank.getLong(o_triggerBits,0); // -->tree

     //if(debug) particleBank.show(); // causes segfault!
     if(debug) { for(int l=0; l<30; l++) printf("-"); printf("\n"); };


     // SEARCH FOR OBSERVABLES, SORTING EACH KIND BY ENERGY
     // ---------------------------------------------------
     //
     // -- reset trajectory-sorting data structures
     for(int p=0; p<nParticles; p++) {
       trajCnt[p] = 0;
       trajArrUS[p]->Clear();
       trajArr[p]->Clear();
       for(int t=0; t<maxTraj; t++) {
         trajSortIdx[p][t] = -1;
         trajE[p][t] = -1;
       };
     };
     
     // -- if it's particle that know about, add its trajectory
     //    to the unsorted trajectory array and its energy to list
     //    of energies for this particle
     particleCntAll = particleList.getSize(); // -->tree
     for(int i=0; i<particleCntAll; i++) {

       // get particle PID and convert it to particle index (in Constants.h)
       pidCur = particleList.getPid(i);
       pIdx = PIDtoIdx(pidCur);

       // skip particles we don't care about (so we don't waste time sorting them)
       if(pIdx==kP || pIdx==kN) pIdx=-10000;

       if(debug) printf(" pid=%d  pIdx=%d\n",pidCur,pIdx);

       if(pIdx>-10000) {

         particleList.getVector4(i,vecObs,PartMass(pIdx));

         // set Trajectory pointer tr to proper allocated Trajectory instance;
         // if there are more instances of this observable than we allocated for,
         // we could allocate more memory for the extra ones (which could cause
         // a memory leak or eventual segfault if maxTraj is not high enough...),
         // or just ignore the extra particles and print a warning 
         if(trajCnt[pIdx]<maxTraj) {
           
           // make sure particle has sensible momentum and energy 
           // (sometimes photons in dnp2018 skim files don't...)
           if( !isnan(vecObs.vect().x()) && !isnan(vecObs.vect().y()) && 
               !isnan(vecObs.vect().z()) && !isnan(vecObs.e()) ) {

             // set tr to allocated Trajectory instance and add it to the unsorted array
             tr = traj[pIdx][trajCnt[pIdx]]; 
             tr->SetMomentum(
               vecObs.vect().x(),
               vecObs.vect().y(),
               vecObs.vect().z()
             );
             trajArrUS[pIdx]->AddLast(tr);

             // add this trajectory's energy to the energy array
             trajE[pIdx][trajCnt[pIdx]] = vecObs.e();

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
     };
     if(debug) printf(">>\n");


     // -- sort each list of each type of observable by E; trajArr will 
     //    be the sorted trajectory array and the 0th entry will be the 
     //    highest energy trajectory 
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
     //
     // ---------------------------------------------------

     
     // PHOTON PAIRING, for pi0s
     // -- this is done by sorting all photon pairs by energy, and picking the
     //    highest-E pair that satisfies basic requirements set in Diphoton class
     //    booleans (e.g., an opening angle cut)

     // if a pi0 has been found in the particle bank, its trajectory will already be in
     // trajArr, and the highest-E pi0 found here will be the one sent into the tree;
     // the user will be warned of this happening
     if(trajCnt[kPi0] > 0) {
       fprintf(stderr,
         "WARNING WARNING WARNING WARNING WARNING: found pi0 in paricle bank!!!\n");
       // note that EventTree is not yet programmed to accept these pions... TODO
     };


     // -- reset data structures for photon pair sorting
     phpCnt = 0;
     for(int php=0; php<phpCntMax; php++) {
       phpIdx[php][0] = -1;
       phpIdx[php][1] = -1;
       phpE[php] = -1;
       phpSortIdx[php] = -1;
     };
     for(int pp=0; pp<maxTraj; pp++) photonUsed[pp] = false;


     // if there are at least 2 photons, we can look for pi0s
     if( trajCnt[kPhoton] >= 2 ) {
       
       if(debugPHPsort) printf(">>>>>>\n");

       // -- loop through pairs of photons
       for(int p0=0; p0<trajCnt[kPhoton]; p0++) {
         for(int p1=p0+1; p1<trajCnt[kPhoton]; p1++) {

           // if there are more than phpCntMax pairs, additional ones are simply ignored
           if(phpCnt < phpCntMax) {

             // fill php index arrays
             phpIdx[phpCnt][0] = p0;
             phpIdx[phpCnt][1] = p1;

             // fill php energy array
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

       // -- sort photon pairs by E
       TMath::Sort(phpCnt,phpE,phpSortIdx);

       // -- test sorting output
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

             // compute kinematics
             diPhotTmp[php]->SetEvent(phot[0],phot[1]);

             // if it satisfies basic requirements, add it to trajArr; if it's the first
             // one added to trajArr (i.e., highest-E), set diphoton tree branches too
             if(diPhotTmp[php]->validDiphoton) {
               for(int h=0; h<2; h++) photonUsed[phpSI[h]] = true; // mark photons used
               trajArr[kDiph]->AddLast(diPhotTmp[php]->Traj);
               if(trajCnt[kDiph]==0) diPhot->SetEvent(phot[0],phot[1]); // -->tree
               trajCnt[kDiph]++;
             };
           };

         };
       };
     };


     

     // analyze hadrons and fill tree
     // ---------------------------------
     
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
   };

   // print benchmark
   reader.showBenchmark();
   printf(" time spend on analysis = %8.5f (ms) events = %8d\n",
     bench.getTime()*1e-9,
     bench.getCounter());

   // write output tree
   tree->Write("tree");

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
   outfile->Close();
};
