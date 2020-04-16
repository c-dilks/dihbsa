// read skim HIPO file and produce a ROOT tree, called an 'outroot file', which contains
// branches associated with dihadrons

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
#include "TRandom.h"


// Clas12Tool
#include "reader.h"
#include "bank.h"
#include "particle.h"
#include "mcparticle.h"
#include "clas12reader.h"

// DihBsa
#include "Constants.h"
#include "DIS.h"
#include "Trajectory.h"
#include "FiducialCuts.h"
#include "Dihadron.h"
#include "Diphoton.h"
#include "EventTree.h"


// methods / vars for MC
Int_t generateHelicity(Int_t idx, Float_t phiH, Float_t phiR); // for MC
Float_t modu(Int_t par, Float_t ph, Float_t pr); // for MC
TRandom * RNG; 
//


int main(int argc, char** argv) {

   // ARGUMENTS
   TString infileN;
   Bool_t augerMode = false;
   if(argc<=1) {
     printf("USAGE: %s [hipo file]\n",argv[0]);
     exit(0);
   };
   if(argc>1) infileN = TString(argv[1]);
   if(argc>2) augerMode = (Bool_t) strtof(argv[2],NULL);
   

   // if true, will only analyze the leading-E hadrons
   // -- should be set to false, so theorists can interpret it
   //    as that is more 'inclusive'
   bool leadingDihOnly = false;


   // debugging flags
   bool debug = 0; // general debugging statements
   bool debugPHPsort = 0; // sorting photon pairs by energy




   // set output ROOT file
   TString outfileN;
   if(augerMode) outfileN = "outroot.root";
   else {
     outfileN = infileN;
     outfileN(TRegexp("^.*/")) = "outroot/";
     outfileN += ".root";
   };
   printf("outfileN = %s\n",outfileN.Data());
   TFile * outfile = new TFile(outfileN,"RECREATE");


   // load libs
   DIS * disEv = new DIS();
   Dihadron * dih = new Dihadron(); dih->useBreit = false;
   Diphoton * diPhot[2];
   for(int h=0; h<2; h++) diPhot[h] = new Diphoton();

   // define trajectories
   Trajectory * had[2]; // for dihadron pair; points to traj instances
   Trajectory * ele; // for electron
   Trajectory * phot[2]; // for diphoton (for pi0)

   // set up structure for sorting particles by E
   //   each particle has a trajectory pointer array,
   //   which will be sorted by energy later
   const Int_t maxTraj = 40;
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
   tree->Branch("eleVertex",disEv->eleVertex,"eleVertex[3]/F");
   tree->Branch("eleStatus",&(disEv->eleStatus),"eleStatus/I");
   tree->Branch("eleChi2pid",&(disEv->eleChi2pid),"eleChi2pid/F");
   
   // electron fiducial cuts branches
   Bool_t eleFidPCAL[FiducialCuts::nLevel];
   Bool_t eleFidDC[FiducialCuts::nLevel];
   TString brsuffix = Form("[%d]",FiducialCuts::nLevel);
   tree->Branch("eleFidPCAL",eleFidPCAL,TString("eleFidPCAL"+brsuffix+"/O"));
   tree->Branch("eleFidDC",eleFidDC,TString("eleFidDC"+brsuffix+"/O"));


   // multiplicity branches
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
   Int_t hadOrder;
   Int_t hadIdx[2]; // particle Idx of each hadron in the pair
   Float_t hadE[2];
   Float_t hadP[2];
   Float_t hadPt[2]; 
   Float_t hadPtq[2];
   Float_t hadEta[2];
   Float_t hadTheta[2];
   Float_t hadPhi[2];
   tree->Branch("pairType",&pairType,"pairType/I"); // pair type number
   tree->Branch("hadOrder",&hadOrder,"hadOrder/I"); // (see hadron pairing algorithm)
   tree->Branch("hadIdx",hadIdx,"hadIdx[2]/I");
   tree->Branch("hadE",hadE,"hadE[2]/F");
   tree->Branch("hadP",hadP,"hadP[2]/F");
   tree->Branch("hadPt",hadPt,"hadPt[2]/F");
   tree->Branch("hadEta",hadEta,"hadEta[2]/F");
   tree->Branch("hadPhi",hadPhi,"hadPhi[2]/F");
   tree->Branch("hadXF",dih->hadXF,"hadXF[2]/F");
   tree->Branch("hadVertex",dih->hadVertex,"hadVertex[2][3]/F");
   tree->Branch("hadStatus",dih->hadStatus,"hadStatus[2]/I");
   tree->Branch("hadChi2pid",dih->hadChi2pid,"hadChi2pid[2]/F");
   /*
   Bool_t hadFidPCAL[2], hadFidDC[2];
   tree->Branch("hadFidPCAL",hadFidPCAL,"hadFidPCAL[2]/O");
   tree->Branch("hadFidDC",hadFidDC,"hadFidDC[2]/O");
   */

   // dihadron branches
   tree->Branch("Mh",&(dih->Mh),"Mh/F");
   tree->Branch("Mmiss",&(dih->Mmiss),"Mmiss/F");
   tree->Branch("Z",dih->z,"Z[2]/F");
   tree->Branch("Zpair",&(dih->zpair),"Zpair/F");
   tree->Branch("xF",&(dih->xF),"xF/F");
   tree->Branch("alpha",&(dih->alpha),"alpha/F");
   tree->Branch("theta",&(dih->theta),"theta/F");
   tree->Branch("zeta",&(dih->zeta),"zeta/F");

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


   // pi0 (diphoton) branches
   // -- they are arrays of length diphCnt; there are three cases:
   //    - diphCnt=0: no diphotons in the dihadron
   //    - diphCnt=1: one dihadron hadrons is a diphoton (pi0 or BG)
   //    - diphCnt=2: both dihadron hadrons are diphotons (pi0-pi0, pi0-BG, or BG-BG)
   Int_t diphCnt,diphCnt_tr;
   Float_t diphPhotE[2][2], diphPhotPt[2][2], diphPhotEta[2][2], diphPhotPhi[2][2];
   Float_t diphPhotVertex[2][2][3];
   Float_t diphPhotChi2pid[2][2];
   Float_t diphE[2], diphZ[2], diphPt[2], diphM[2], diphAlpha[2], diphEta[2], diphPhi[2];
   tree->Branch("diphCnt",&diphCnt_tr,"diphCnt/I"); // number of diphotons {0,1,2}
   tree->Branch("diphPhotE",diphPhotE,"diphPhotE[diphCnt][2]/F"); // photon energy
   tree->Branch("diphPhotPt",diphPhotPt,"diphPhotPt[diphCnt][2]/F"); // photon pT
   tree->Branch("diphPhotEta",diphPhotEta,"diphPhotEta[diphCnt][2]/F"); // photon eta
   tree->Branch("diphPhotPhi",diphPhotPhi,"diphPhotPhi[diphCnt][2]/F"); // photon phi
   tree->Branch("diphPhotVertex",diphPhotVertex,"diphPhotVertex[diphCnt][2][3]/F"); // photon vertex
   tree->Branch("diphPhotChi2pid",diphPhotChi2pid,"diphPhotChi2pid[diphCnt][2]/F"); // photon pid chi2
   tree->Branch("diphE",diphE,"diphE[diphCnt]/F"); // energy
   tree->Branch("diphZ",diphZ,"diphZ[diphCnt]/F"); // energy imbalance
   tree->Branch("diphPt",diphPt,"diphPt[diphCnt]/F"); // pT
   tree->Branch("diphM",diphM,"diphM[diphCnt]/F"); // invariant mass
   tree->Branch("diphAlpha",diphAlpha,"diphAlpha[diphCnt]/F"); // opening angle
   tree->Branch("diphEta",diphEta,"diphEta[diphCnt]/F"); // eta
   tree->Branch("diphPhi",diphPhi,"diphPhi[diphCnt]/F"); // pT

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
   clas12::clas12reader reader(infileN.Data());

   hipo::dictionary factory;
   reader.getReader().readDictionary(factory);
   //factory.show();
   hipo::event readerEvent;
   hipo::bank mcParticle(factory.getSchema("MC::Particle"));


   // define observable variables
   TLorentzVector vecObs; // observable 4-momentum
   Float_t vecObsP[3];
   Float_t vertex[3];
   Float_t chi2pid;
   Int_t status;


   Int_t pidCur,pIdx;
   Int_t i1,i2;

   Bool_t pairsExist;
   Bool_t printPi0Warning = true;


   // MC: variables
   Bool_t printWarning = true;
   RNG = new TRandom(14972);
   TString genfileN;
   Bool_t genSuccess,foundMatch;
   Int_t genPID;
   Float_t hadMD[2];
   Float_t hadMDmin[2];
   Float_t genEleP, genElePmax;
   TLorentzVector genHadVec[2];
   TLorentzVector genHadVecTmp;
   TLorentzVector genEleVecTmp, genEleVec;
   Trajectory * genEleTraj = new Trajectory(kE);
   Trajectory * genHadTraj[2];
   for(int h=0; h<2; h++) genHadTraj[h] = new Trajectory(kPhoton); // (Idx is unimportant)
   DIS * genDIS = new DIS();
   Dihadron * genDih = new Dihadron();


   // determine MC mode, using PARTICLE_BANK macro, defined in config.mk
   Bool_t MCgenMode, MCrecMode;
   Int_t whichEle; // set to DIS electron (highest-E electron, but if reading MC::Lund,
                   // it's the 2nd highest-E)
#if PARTICLE_BANK == 0 // REC::Particle -- for data
   MCgenMode = false;
   MCrecMode = false;
   whichEle = 0;
#elif PARTICLE_BANK == 1 // MC::Lund -- for MC generated
   MCgenMode = true;
   MCrecMode = false;
   whichEle = 1;
#elif PARTICLE_BANK == 2 // MC::Particle -- for MC generated
   MCgenMode = true;
   MCrecMode = false;
   whichEle = 0;
#elif PARTICLE_BANK == 3 // REC::Particle -- for MC reconstructed
   MCgenMode = false;
   MCrecMode = true;
   whichEle = 0;
#else
  fprintf(stderr,"ERROR: preprocessor macro PARTICLE_BANK=%d is undefined\n",
    PARTICLE_BANK);
  exit(0);
#endif


   // MC branches
   Int_t helicityMC[EventTree::NhelicityMC];
   TString brStr;
   Float_t MD;
   Float_t gen_eleE, gen_elePt, gen_eleEta, gen_elePhi;
   Float_t gen_hadE[2];
   Float_t gen_hadPt[2];
   Float_t gen_hadEta[2];
   Float_t gen_hadPhi[2];
   if(MCrecMode || MCgenMode) {
     brStr = Form("helicityMC[%d]/I",EventTree::NhelicityMC);
     tree->Branch("helicityMC",helicityMC,brStr);
   };
   if(MCrecMode) {
     tree->Branch("matchDiff",&MD,"matchDiff/F");
     tree->Branch("gen_eleE",&gen_eleE,"gen_eleE/F");
     tree->Branch("gen_elePt",&gen_elePt,"gen_elePt/F");
     tree->Branch("gen_eleEta",&gen_eleEta,"gen_eleEta/F");
     tree->Branch("gen_elePhi",&gen_elePhi,"gen_elePhi/F");
     tree->Branch("gen_hadE",gen_hadE,"gen_hadE[2]/F");
     tree->Branch("gen_hadPt",gen_hadPt,"gen_hadPt[2]/F");
     tree->Branch("gen_hadEta",gen_hadEta,"gen_hadEta[2]/F");
     tree->Branch("gen_hadPhi",gen_hadPhi,"gen_hadPhi[2]/F");
   };



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
   Int_t itCount = 0;

   int end1,end2;


   // ----------------------------------------------------
   // EVENT LOOP
   // ----------------------------------------------------
   printf("begin event loop...\n");

   while(reader.next()==true) {

     //if(evCount>100) { fprintf(stderr,"BREAKING LOOP HERE!!!\n\n"); break; };
     //itCount++; if(itCount>10) { fprintf(stderr,"BREAKING LOOP HERE!!!\n\n"); break; };

     bench.resume();


     // reset branches
     disEv->ResetVars();
     dih->ResetVars();
     for(int h=0; h<2; h++) {
       hadIdx[h] = UNDEF;
       hadE[h] = UNDEF;
       hadP[h] = UNDEF;
       hadPt[h] = UNDEF;
       hadPtq[h] = UNDEF;
       hadEta[h] = UNDEF;
       hadTheta[h] = UNDEF;
       hadPhi[h] = UNDEF;
     };
     pairType = UNDEF;
     for(int h=0; h<2; h++) diPhot[h]->ResetVars();
     for(int php=0; php<phpCntMax; php++) diPhotTmp[php]->ResetVars();
     diphCnt = 0;



     // read event-header stuff
     evnum = reader.runconfig()->getEvent(); // -->tree
     runnum = reader.runconfig()->getRun(); // -->tree
     triggerBits = reader.runconfig()->getTrigger(); // -->tree
     torus = reader.runconfig()->getTorus(); // -->tree
     solenoid = reader.runconfig()->getSolenoid(); // -->tree
     //printf("schema name = %s\n",reader.head()->getSchema().getName().data());


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


     // HELICITY
     helicity = reader.event()->getHelicity(); // -->tree



     
     // ---------------------------------------------------
     // PARTICLE LOOP
     // ---------------------------------------------------
     // -- read in each particle and put them into trajArr, which will be sorted
     //    afterward
     

#if PARTICLE_BANK == 0 // REC::Particle for reconstructed particles from data
     particleCntAll = reader.getNParticles(); // -->tree
     if(debug) printf("reader.getNParticles() = %d\n",particleCntAll);
     for(auto & part : reader.getDetParticles()) {
       pidCur = part->getPid();
       vecObsP[eX] = part->par()->getPx();
       vecObsP[eY] = part->par()->getPy();
       vecObsP[eZ] = part->par()->getPz();
       vertex[eX] = part->par()->getVx();
       vertex[eY] = part->par()->getVy();
       vertex[eZ] = part->par()->getVz();
       chi2pid = part->par()->getChi2Pid();
       status = part->par()->getStatus();
#elif PARTICLE_BANK == 1 // MC::Lund for reading MC-generated particles
     particleCntAll = (reader.mcparts())->getRows(); // -->tree
     for(int rr=0; rr<particleCntAll; rr++) {
       pidCur = (reader.mcparts())->getPid(rr);
       vecObsP[eX] = (reader.mcparts())->getPx(rr);
       vecObsP[eY] = (reader.mcparts())->getPy(rr);
       vecObsP[eZ] = (reader.mcparts())->getPz(rr);
       vertex[eX] = (reader.mcparts())->getVx(rr);
       vertex[eY] = (reader.mcparts())->getVy(rr);
       vertex[eZ] = (reader.mcparts())->getVz(rr);
       chi2pid = UNDEF;
       status = 0;
#elif PARTICLE_BANK == 2 // MC::Particle for reading MC-generated particles
     reader.getReader().read(readerEvent);
     readerEvent.getStructure(mcParticle);
     particleCntAll = mcParticle.getRows(); // -->tree
     for(int rr=0; rr<particleCntAll; rr++) {
       pidCur = mcParticle.getInt("pid",rr);
       vecObsP[eX] = mcParticle.getFloat("px",rr);
       vecObsP[eY] = mcParticle.getFloat("py",rr);
       vecObsP[eZ] = mcParticle.getFloat("pz",rr);
       vertex[eX] = mcParticle.getFloat("vx",rr);
       vertex[eY] = mcParticle.getFloat("vy",rr);
       vertex[eZ] = mcParticle.getFloat("vz",rr);
       chi2pid = UNDEF;
       status = 0;
#elif PARTICLE_BANK == 3 // REC::Particle for reconstructed particles from MC
     reader.getReader().read(readerEvent);
     readerEvent.getStructure(mcParticle); // for matching to MCgen
     particleCntAll = reader.getNParticles(); // -->tree
     if(debug) printf("reader.getNParticles() = %d\n",particleCntAll);
     for(auto & part : reader.getDetParticles()) {
       pidCur = part->getPid();
       vecObsP[eX] = part->par()->getPx();
       vecObsP[eY] = part->par()->getPy();
       vecObsP[eZ] = part->par()->getPz();
       vertex[eX] = part->par()->getVx();
       vertex[eY] = part->par()->getVy();
       vertex[eZ] = part->par()->getVz();
       chi2pid = part->par()->getChi2Pid();
       status = part->par()->getStatus();
#endif

       

       // convert PID to local particle index; if it's not defined in Constants.h, pIdx
       // will be UNDEF and this particle will be ignored
       pIdx = PIDtoIdx(pidCur);
       if(pIdx==kP || pIdx==kN) pIdx=UNDEF; // also skip protons and neutrons
       

       // if it's an electron with chi2pid==0.00000, it may not be *the* scattered
       // electron, so we skip it
       if(pIdx==kE && fabs(chi2pid)<1e-5) pIdx=UNDEF;


       if(pIdx>UNDEF) {
         if(debug) printf("\nNEXT PARTICLE: --> pid=%d  pIdx=%d\n",pidCur,pIdx);
        
         // set Trajectory pointer tr to proper allocated Trajectory instance; if there
         // are more instances of this observable than we allocated for, (this max
         // number is "maxTraj"), just ignore the additional particles and print a
         // warning to stderr
         if(trajCnt[pIdx]<maxTraj) {
           
           // make sure particle has sensible momentum (!=NaN)
           // (sometimes photons in dnp2018 skim files don't...)
           if( !isnan(vecObsP[eX]) && !isnan(vecObsP[eY]) && !isnan(vecObsP[eZ]) ) {

             // set tr to allocated Trajectory instance
             tr = traj[pIdx][trajCnt[pIdx]]; 

             // set tr kinematics
             vecObs.SetXYZM(vecObsP[eX],vecObsP[eY],vecObsP[eZ],PartMass(pIdx));
             tr->SetVec(vecObs);
             tr->SetVertex(vertex[eX],vertex[eY],vertex[eZ]);
             tr->chi2pid = chi2pid;
             tr->Status = status;
             
             // set tr FiducialCuts info (note: Trajectory derives from FiducialCuts)
#if PARTICLE_BANK == 0 || PARTICLE_BANK == 3
             tr->enableFiducialCut = true; // (default enableFiducialCut is false)
             tr->torus = torus;
             // -- PCAL from REC::Calorimeter
             tr->pcalSec = part->cal(clas12::PCAL)->getSector();
             tr->pcalLayer = part->cal(clas12::PCAL)->getLayer();
             tr->pcalL[FiducialCuts::u] = part->cal(clas12::PCAL)->getLu();
             tr->pcalL[FiducialCuts::v] = part->cal(clas12::PCAL)->getLv();
             tr->pcalL[FiducialCuts::w] = part->cal(clas12::PCAL)->getLw();
             // note: needed to add getLw implmentation in Clas12Tool's
             //       clas12::calorimeter, it was missing in my current version
             // -- DC from REC::Track
             tr->dcTrackDetector = part->trk(clas12::DC)->getDetector();
             tr->dcSec = part->trk(clas12::DC)->getSector();
             // -- DC from REC::Traj
             //    'r' loops over three regions (layers 6,18,30)
             //    regLayer({0,1,2}) returns {6,18,30}
             for(int r=0; r<FiducialCuts::nReg; r++) {
               tr->dcTrajDetector[r] = 
                 part->traj(clas12::DC,FiducialCuts::regLayer(r))->getDetector();
               tr->dcTrajLayer[r] = 
                 part->traj(clas12::DC,FiducialCuts::regLayer(r))->getLayer();
               tr->dcTraj[r][FiducialCuts::x] = 
                 part->traj(clas12::DC,FiducialCuts::regLayer(r))->getX();
               tr->dcTraj[r][FiducialCuts::y] = 
                 part->traj(clas12::DC,FiducialCuts::regLayer(r))->getY();
               tr->dcTraj[r][FiducialCuts::z] = 
                 part->traj(clas12::DC,FiducialCuts::regLayer(r))->getZ();
             };
             // -- print
             if(debug) tr->PrintFiducialCuts(FiducialCuts::cutMedium);
#endif

             // add tr to unsorted Trajectory array, and energy to the energy array
             trajArrUS[pIdx]->AddLast(tr);
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
     if(debug) Tools::PrintSeparator(30,">");
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
       if(printPi0Warning) {
         fprintf(stderr,
           "WARNING WARNING WARNING WARNING WARNING: found pi0 in paricle bank!!!\n");
         printPi0Warning = false;
       };
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
               //
               // aqui - TODO - this will need to be made more inclusive
               //
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
           diphPhotChi2pid[dp][h] = diPhot[dp]->photChi2pid[h]; // -->tree
           for(int c=0; c<3; c++) {
             diphPhotVertex[dp][h][c] = diPhot[dp]->photVertex[h][c]; // -->tree
           };
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

     if(debug) {
       printf("sorted list: %s cnt=%d\n",PartName(kDiph).Data(),trajCnt[kDiph]);
       for(int t=0; t<trajCnt[kDiph]; t++)
         ((Trajectory*)(trajArr[kDiph]->At(t)))->Vec.Print();
     };

     

     // ---------------------------------------------------
     // HADRON PAIRING, and fill the tree
     // ---------------------------------------------------
     

     // first make sure there's a scattered electron
     if(trajCnt[kE] > whichEle) {

       // compute DIS kinematics
       ele = (Trajectory*) trajArr[kE]->At(whichEle); 
       disEv->SetElectron(ele);
       disEv->Analyse(); // -->tree


       // evaluate fiducial cuts for electron
       // (force them to be true for MC-generated particles)
#if PARTICLE_BANK == 0 || PARTICLE_BANK == 3 // REC::Particle (data or MC)
       for(int l=0; l<FiducialCuts::nLevel; l++) {
         eleFidPCAL[l] = ele->FidPCAL(l); // -->tree
         eleFidDC[l] = ele->FidDC(l); // -->tree
       };
#elif PARTICLE_BANK == 1 || PARTICLE_BANK == 2 // MC generated banks: just set to true
       for(int l=0; l<FiducialCuts::nLevel; l++) {
         eleFidPCAL[l] = true; // -->tree
         eleFidDC[l] = true; // -->tree
       };
#endif


       // look for "observable pairs" -- these are pairs that are used to form
       // form dihadrons; only the desired observables are paired (see observable_enum
       // in src/Constants.h)
       for(Int_t o1=0; o1<nObservables; o1++) {
         for(Int_t o2=0; o2<nObservables; o2++) {

           // each dihadron channel has a specific order, e.g., pi+pi- and pi-pi+ are
           // the same channel, but the order pi+pi- is preferred. There are two tools
           // in src/Constants.h that help acheive this: dihHadIdx and IterPair
           // - dihHadIdx takes two particle indexes, along with 'h' = qA or qB, and 
           //   returns the particle index that should be associated with 'h'
           //   -> example: dihHadIdx(kPip,kPim,qA) returns kPip
           //               dihHadIdx(kPip,kPim,qB) returns kPim
           //          also dihHadIdx(kPim,kPip,qA) returns kPip
           // - IterPair takes a pair of observable indices (here o1 and o2), and
           //   returns 'true' if this pair is in the proper order; it also converts
           //   o1 and o2 to particle indices 
           //   -> example: 
           //      IterPair(sPip,sPim,i1,i2) returns true,  with i1=kPip and i2=kPim
           //      IterPair(sPim,sPip,i1,i2) returns false, with i1=kPim and i2=kPip
           //   -> the return value of IterPair is useful within this loop over o1 and
           //      o2 to ensure we only look at unique pairs (and don't double count pi+pi-
           //      with pi-pi+)


           // convert observable index to particle index;
           // IterPair returns true if o1 and o2 follow conventional dihadron ordering
           // (see above for wall of text for details)
           if(IterPair(o1,o2,i1,i2)) {

             // now i1 is set to Constants::OI(o1), and likewise for i2


             for(int h=0; h<2; h++) hadIdx[h] = dihHadIdx(i1,i2,h); // -->tree
             pairType = EncodePairType(hadIdx[qA],hadIdx[qB]); // -->tree



             // HADRON PAIRING ALGORITHM
             // ---------------------------------------
             // for this dihadron channel, denoted (p,q), we must pair up all the p's
             // with all the q's; this is acheived by another nested loop

             // define loop ends:
             // - consider lists of hadrons (p1,p2,p3) and (q1,q2,q3) (ordered by energy)
             // - if the hadron species are different, loop through each list completely
             //   so that every possible pair is considered
             //   -> the above lists pair as p1q1, p1q2, p1q3
             //                              p2q1, p2q2, p2q3
             //                              p3q1, p3q2, p3q3
             //   -> if leadingDihOnly==true, only take p1q1
             //   -> outer loop and inner loop both begin at 0
             //   -> outer loop and inner loop both end at trajCnt
             // - if the hadron species are the same, loop such that pairs aren't
             //   double-counted; the folowing matrix is the same as above, but
             //   with reduntant pairs replaced by xxxx:
             //   ->                         xxxx, p1p2, p1p3
             //                              xxxx, xxxx, p2p3
             //                              xxxx, xxxx, xxxx
             //   -> if leadingDihOnly==true, only take p1p2
             //   -> outer loop begins at 0;
             //      inner loop begins at outer loop iterator+1 (first column is all xxxx)
             //   -> outer loop ends at trajCnt-1 (last row is all xxxx),
             //      inner loop ends at trajCnt
             //   
             // - hadOrder: this is a number that iterates for each hadron pair
             //   - if hadOrder==0, you can safely assume this is the leading hadron pair
             //   - for hadOrder>0, it is just an iterator for the nested outer/inner
             //     loops, so the remaining dihadrons are not ordered by energy, but
             //     instead by how they are paired
             //   
             //
             if(i1!=i2) {
               end1 = leadingDihOnly ? 1 : trajCnt[i1];
               end2 = leadingDihOnly ? 1 : trajCnt[i2];
               pairsExist = trajCnt[i1]>=1 && trajCnt[i2]>=1;
             }
             else {
               end1 = leadingDihOnly ? 1 : trajCnt[i1] - 1;
               end2 = leadingDihOnly ? 2 : trajCnt[i2];
               pairsExist = trajCnt[i1]>=2;
             };

             hadOrder = 0;


             if(pairsExist) {
               // outer loop, for hadron qA (0)
               for(int hh1=0; hh1<end1; hh1++) {
                 // inner loop, for hadron qB (1)
                 for(int hh2=( i1==i2 ? hh1+1 : 0 ); hh2<end2; hh2++) {

                   if(debug) {
                     printf("\n>>> BEGIN DIHADRON EVENT %s\n",PairName(i1,i2).Data());
                     printf("    hadOrder = %d\n",hadOrder);
                   };



                   // get hadron kinematics
                   had[qA] = (Trajectory*) trajArr[hadIdx[qA]]->At(hh1);
                   had[qB] = (Trajectory*) trajArr[hadIdx[qB]]->At(hh2);

                   for(int h=0; h<2; h++) {
                     hadE[h] = (had[h]->Vec).E(); // -->tree
                     hadP[h] = (had[h]->Vec).P(); // -->tree
                     hadPt[h] = (had[h]->Vec).Pt(); // -->tree
                     hadPtq[h] = had[h]->Ptq(disEv->vecQ); // -->tree

                     if(hadE[h]>0 && hadPt[h]>0) {
                       hadEta[h] = (had[h]->Vec).Eta(); // -->tree
                       hadTheta[h] = (had[h]->Vec).Theta();
                       hadPhi[h] = (had[h]->Vec).Phi(); // -->tree
                     } else {
                       hadEta[h] = UNDEF;
                       hadTheta[h] = UNDEF;
                       hadPhi[h] = UNDEF;
                     };

                     /*
                     hadFidPCAL[h] = had[h]->FidPCAL(FiducialCuts::cutMedium); // -->tree
                     hadFidDC[h] = had[h]->FidDC(FiducialCuts::cutMedium); // -->tree
                     */

                     if(debug) {
                       printf("[+] %s 4-momentum:\n",(had[h]->Title()).Data());
                       (had[h]->Vec).Print();
                     };
                   };

                   /*
                   dih->debugTheta = pairType==0x34 ? true : false;
                   if(dih->debugTheta) printf("\nevnum = %d\n",evnum);
                   */

                   // compute dihadron kinematics
                   dih->SetEvent(had[qA],had[qB],disEv); // -->tree

                   // set diphCnt_tr to corrected diphCnt, for filling diphoton branches
                   // -- aqui - TODO - fix this too!
                   // -- this solves a subtle case: assume we have a pi+ and two
                   // diphotons; this will resolve to two dihadrons: A=(pi+,diphoton) and
                   // B=(diphoton,diphoton), but diphCnt=2 for this event for both. We
                   // thus need to force diphCnt to 1 for dihadron A, and to 2 for
                   // dihadron B
                   diphCnt_tr = 0; // -->tree
                   if(diphCnt>0) {
                     if(hadIdx[qA]==kDiph || hadIdx[qB]==kDiph) diphCnt_tr = 1; // -->tree
                     if(hadIdx[qA]==kDiph && hadIdx[qB]==kDiph) diphCnt_tr = 2; // -->tree
                   };



                   // MC helicity injection
                   if(MCgenMode || MCrecMode) {

                     if(MCgenMode && MCrecMode) {
                       fprintf(stderr,"ERROR: both MCgenMode and MCrecMode are true\n");
                       return 0;
                     };


                     if(MCgenMode) {
                       // generate helicity based on injected asymmetry
                       for(int hh=0; hh<EventTree::NhelicityMC; hh++) {
                         helicityMC[hh] = generateHelicity(hh,dih->PhiH,dih->PhiR);
                       };
                     };

                     // use matching MCgen kinematics in helicity generation
                     /*
                     if(MCrecMode) {

                       // reset variables used for generated hadron matching
                       for(int hp=0; hp<2; hp++) {
                         genHadVec[hp].SetXYZM(0.0,0.0,0.0,0.0);
                         genDIS->ResetVars();
                         genDih->ResetVars();
                         hadMDmin[hp] = 1e6;
                       };
                       genEleVec.SetXYZM(0.0,0.0,0.0,0.0);
                       genElePmax = 0;

                       // loop through MC::Particle and find closest matching hadrons
                       // as well as the generated scattered electron
                       for(int mp=0; mp<mcParticle.getRows(); mp++) {
                         genPID = mcParticle.getInt("pid",mp);
                         
                         // generated scattered electron
                         if(genPID==PartPID(kE)) {
                           genEleVecTmp.SetXYZM(
                             mcParticle.getFloat("px",mp),
                             mcParticle.getFloat("py",mp),
                             mcParticle.getFloat("pz",mp),
                             PartMass(kE)
                           );
                           genEleP = genEleVecTmp.P();
                           if(genEleP > genElePmax) {
                             genElePmax = genEleP;
                             genEleVec.SetXYZM(
                               genEleVecTmp.Px(),
                               genEleVecTmp.Py(),
                               genEleVecTmp.Pz(),
                               genEleVecTmp.M()
                             );
                           };
                         }
                         
                         // generated hadrons
                         else {
                           for(int hp=0; hp<2; hp++) {
                             if(genPID==PartPID(hadIdx[hp])) {

                               genHadVecTmp.SetXYZM(
                                 mcParticle.getFloat("px",mp),
                                 mcParticle.getFloat("py",mp),
                                 mcParticle.getFloat("pz",mp),
                                 PartMass(hadIdx[hp])
                               );

                               // calculate 2-dim Euclidean distance ("matchDiff", or
                               // "MD") between MCrec and MCgen hadrons; the smallest MD
                               // is taken as the "best matched" hadron; we can cut on
                               // MD later on in order to select only the matches we
                               // have confidence in
                               hadMD[hp] = TMath::Sqrt(
                                 TMath::Power(
                                   Tools::AdjAngle(hadTheta[hp] - genHadVecTmp.Theta()),
                                   2) +
                                 TMath::Power(
                                   Tools::AdjAngle(hadPhi[hp] - genHadVecTmp.Phi()),
                                   2)
                               );

                               if(hadMD[hp] < hadMDmin[hp]) {
                                 hadMDmin[hp] = hadMD[hp];
                                 genHadVec[hp].SetXYZM(
                                   genHadVecTmp.Px(),
                                   genHadVecTmp.Py(),
                                   genHadVecTmp.Pz(),
                                   genHadVecTmp.M()
                                 );
                               };

                             };
                           };
                         };
                       };
                       
                       // if we found a generated scattered electron and dihadron,
                       // compute kinematics and assign helicities
                       if(genElePmax>0 && hadMDmin[qA]<10 && hadMDmin[qB]<10) {
                         genEleTraj->SetVec(genEleVec);
                         for(int hp=0; hp<2; hp++) {
                           genHadTraj[hp]->SetVec(genHadVec[hp]);
                           genHadTraj[hp]->SetIdx(hadIdx[hp]);
                         };
                         genDIS->SetElectron(genEleTraj);
                         genDIS->Analyse();
                         genDih->SetEvent(genHadTraj[qA],genHadTraj[qB],genDIS);

                         // generate helicity based on injected asymmetry
                         // - this uses the generated dihadron kinematics
                         for(int hh=0; hh<EventTree::NhelicityMC; hh++) {
                           helicityMC[hh] = generateHelicity(
                             hh, genDih->PhiH, genDih->PhiR );
                         };

                         // set MCgenMatched branches
                         MD = TMath::Sqrt(
                           TMath::Power(hadMDmin[qA],2) +
                           TMath::Power(hadMDmin[qB],2) ); // 4-dim. matchDiff
                         gen_eleE = (genEleTraj->Vec).E();
                         gen_elePt = (genEleTraj->Vec).Pt();
                         gen_eleEta = (genEleTraj->Vec).Eta();
                         gen_elePhi = (genEleTraj->Vec).Phi();
                         for(int hp=0; hp<2; hp++) {
                           gen_hadE[hp] = (genHadVec[hp]).E();
                           gen_hadPt[hp] = (genHadVec[hp]).Pt();
                           gen_hadEta[hp] = (genHadVec[hp]).Eta();
                           gen_hadPhi[hp] = (genHadVec[hp]).Phi();
                         };
                       } else {
                         // if no match candidate was found, set everything to undefined
                         for(int hh=0; hh<EventTree::NhelicityMC; hh++) helicityMC[hh]=0;
                         MD = UNDEF;
                         gen_eleE = UNDEF;
                         gen_elePt = UNDEF;
                         gen_eleEta = UNDEF;
                         gen_elePhi = UNDEF;
                         for(int hp=0; hp<2; hp++) {
                           gen_hadE[hp] = UNDEF;
                           gen_hadPt[hp] = UNDEF;
                           gen_hadEta[hp] = UNDEF;
                           gen_hadPhi[hp] = UNDEF;
                         };
                       };

                     };
                     */


                     // use MCrec kinematics in helicity generation (no matching)
                     ///*
                     if(MCrecMode) {
                       for(int hh=0; hh<EventTree::NhelicityMC; hh++) {
                         helicityMC[hh] = generateHelicity(hh,dih->PhiH,dih->PhiR);
                       };
                       // set MCgen kinematics to MCrec kinematics and MD to 0, so that
                       // all MCrec dihadrons pass matching cuts
                       MD = 0;
                       gen_eleE = disEv->eleE;
                       gen_elePt = disEv->elePt;
                       gen_eleEta = disEv->eleEta;
                       gen_elePhi = disEv->elePhi;
                       for(int hp=0; hp<2; hp++) {
                         gen_hadE[hp] = hadE[hp];
                         gen_hadPt[hp] = hadPt[hp];
                         gen_hadEta[hp] = hadEta[hp];
                         gen_hadPhi[hp] = hadPhi[hp];
                       };
                     };
                     //*/


                     helicity = 0; // set data helicity to zero
                     // ... and print a warning 
                     // (to prevent from accidentally altering the helicity of real data)
                     if(printWarning) {
                       fprintf(stderr,"WARNING WARNING WARNING: ");
                       fprintf(stderr,"helicity has been altered (likely for MC study)!!!!!\n");
                       printWarning = false;
                     };
                   };

                   // fill tree
                   // aqui - TODO - if there's a diphoton, the tree will not fill unless
                   // it's the leading-E diphoton (because at the moment we cannot trust
                   // the kinematics for subleading-E diphotons)
                   if( diphCnt_tr==0 || (diphCnt_tr==1 && hadOrder==0) ) {
                     tree->Fill();
                     if(debug) printf("write to tree\n");
                   };

                   // increment event counter
                   evCount++;
                   if(evCount%100000==0) printf("[---] %d dihadron events found\n",evCount);

                   // increment hadOrder iterator
                   hadOrder++; //-->tree

                 }; // eo hadron pairing inner loop (iterator hh2)
               }; // eo hadron pairing outer loop (iterator hh1)
             }; // eo if pairsExist
             if(debug && !pairsExist) 
               printf("no %s dihadrons.\n",PairName(i1,i2).Data());
           }; // eo if(IterPair)
         }; // eo for o2
       }; // eo for o1
       if(debug) printf("\n");
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
   outfile->cd();
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
   outfile->Close();
};


Int_t generateHelicity(Int_t idx, Float_t phiH, Float_t phiR) {

  Float_t amp = 0.10;
  Float_t asym;

  switch(idx) {
    case 0: asym = 0; break;

    case 1: asym = amp*modu(1,phiH,phiR); break;
    case 2: asym = amp*modu(2,phiH,phiR); break;
    case 3: asym = amp*modu(3,phiH,phiR); break;
    case 4: asym = amp*modu(4,phiH,phiR); break;

    case 5: asym = amp*modu(1,phiH,phiR) + (amp+0.02)*modu(3,phiH,phiR); break;
    case 6: asym = amp*modu(2,phiH,phiR) + (amp+0.02)*modu(3,phiH,phiR); break;
    case 7: asym = amp*modu(3,phiH,phiR) + (amp+0.02)*modu(3,phiH,phiR); break;
    case 8: asym = amp*modu(4,phiH,phiR) + (amp+0.02)*modu(3,phiH,phiR); break;

    case 9:  asym = amp*modu(1,phiH,phiR) - amp*modu(2,phiH,phiR); break;
    case 10: asym = amp*modu(1,phiH,phiR) - amp*modu(2,phiH,phiR) + (amp+0.02)*modu(3,phiH,phiR); break;
    case 11: asym = 0.2*modu(1,phiH,phiR) + 0.2*modu(3,phiH,phiR) + 0.1*modu(2,phiH,phiR); break; // (to compare with Timothy)
    /* EventTree::NhelicityMC must equal the number of cases */
    default: fprintf(stderr,"generateHelicity undefined for idx=%d\n",idx); return 0;
  };

  asym *= 0.86; // polarization factor

  // generate random number in [0,1]
  Float_t rn = RNG->Uniform();

  // return helicity assignment:  2 = spin -   3 = spin +
  return ( rn < 0.5*(1+asym) )  ?  3 : 2;
};
Float_t modu(Int_t par, Float_t ph, Float_t pr) {
  switch(par) {
    case 1: return TMath::Sin(pr); break;
    case 2: return TMath::Sin(ph-pr); break;
    case 3: return TMath::Sin(ph); break;
    case 4: return TMath::Sin(2*ph-pr); break;
    default: fprintf(stderr,"modu unknown\n"); return 0;
  };
};




