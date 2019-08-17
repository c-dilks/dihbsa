#include "EventTree.h"

ClassImp(EventTree)

using namespace std;


EventTree::EventTree(TString filelist, Int_t whichPair_) {
  printf("EventTree instantiated\n");

  debug = true;

  useDiphBG = false;

  DecodePairType(whichPair_,whichHad[qA],whichHad[qB]);
  printf("\n>>> DIHADRON SELECTION: %s\n\n",PairName(whichHad[qA],whichHad[qB]).Data());


  printf("reading tree chain from %s\n",filelist.Data());
  chain = new TChain("tree");
  chain->Add(filelist);

  ENT = chain->GetEntries();
  printf("number of entries: %lld\n",ENT);

  chain->SetBranchAddress("W",&W);
  chain->SetBranchAddress("Q2",&Q2);
  chain->SetBranchAddress("Nu",&Nu);
  chain->SetBranchAddress("x",&x);
  chain->SetBranchAddress("y",&y);

  chain->SetBranchAddress("eleE",&eleE);
  chain->SetBranchAddress("eleP",&eleP);
  chain->SetBranchAddress("elePt",&elePt);
  chain->SetBranchAddress("eleEta",&eleEta);
  chain->SetBranchAddress("elePhi",&elePhi);

  chain->SetBranchAddress("pairType",&pairType);
  chain->SetBranchAddress("hadIdx",hadIdx);
  chain->SetBranchAddress("hadE",hadE);
  chain->SetBranchAddress("hadP",hadP);
  chain->SetBranchAddress("hadPt",hadPt);
  chain->SetBranchAddress("hadEta",hadEta);
  chain->SetBranchAddress("hadPhi",hadPhi);

  // (these branches, which were added for cross-checking, 
  //  aren't yet in current ROOT files)
  if(chain->GetBranch("hadPtq")) chain->SetBranchAddress("hadPtq",hadPtq);
  else { for(int h=0; h<2; h++) hadPtq[h]=-10000; };
  if(chain->GetBranch("hadXF")) chain->SetBranchAddress("hadXF",hadXF);
  else { for(int h=0; h<2; h++) hadXF[h]=-10000; };
  /////////////

  chain->SetBranchAddress("particleCnt",particleCnt);
  chain->SetBranchAddress("particleCntAll",&particleCntAll);

  chain->SetBranchAddress("Mh",&Mh);
  chain->SetBranchAddress("Mmiss",&Mmiss);
  chain->SetBranchAddress("Z",Z);
  chain->SetBranchAddress("Zpair",&Zpair);
  chain->SetBranchAddress("xF",&xF);
  chain->SetBranchAddress("alpha",&alpha);
  chain->SetBranchAddress("theta",&theta);

  chain->SetBranchAddress("Ph",&Ph);
  chain->SetBranchAddress("PhPerp",&PhPerp);
  chain->SetBranchAddress("PhEta",&PhEta);
  chain->SetBranchAddress("PhPhi",&PhPhi);
  chain->SetBranchAddress("R",&R);
  chain->SetBranchAddress("RPerp",&RPerp);
  chain->SetBranchAddress("RT",&RT);

  chain->SetBranchAddress("PhiH",&PhiH);

  chain->SetBranchAddress("PhiRq",&PhiRq);
  chain->SetBranchAddress("PhiRp",&PhiRp);
  chain->SetBranchAddress("PhiRp_r",&PhiRp_r);
  chain->SetBranchAddress("PhiRp_g",&PhiRp_g);

  /*
  chain->SetBranchAddress("b_PhiRq",&b_PhiRq);
  chain->SetBranchAddress("b_PhiRp",&b_PhiRp);
  chain->SetBranchAddress("b_PhiRp_r",&b_PhiRp_r);
  chain->SetBranchAddress("b_PhiRp_g",&b_PhiRp_g);
  */

  chain->SetBranchAddress("runnum",&runnum);
  chain->SetBranchAddress("evnum",&evnum);
  chain->SetBranchAddress("helicity",&helicity);
  chain->SetBranchAddress("torus",&torus);
  chain->SetBranchAddress("triggerBits",&triggerBits);

  chain->SetBranchAddress("diphCnt",&diphCnt);
  chain->SetBranchAddress("diphPhotE",diphPhotE);
  chain->SetBranchAddress("diphPhotPt",diphPhotPt);
  chain->SetBranchAddress("diphPhotEta",diphPhotEta);
  chain->SetBranchAddress("diphPhotPhi",diphPhotPhi);
  chain->SetBranchAddress("diphE",diphE);
  chain->SetBranchAddress("diphZ",diphZ);
  chain->SetBranchAddress("diphPt",diphPt);
  chain->SetBranchAddress("diphM",diphM);
  chain->SetBranchAddress("diphAlpha",diphAlpha);
  chain->SetBranchAddress("diphEta",diphEta);
  chain->SetBranchAddress("diphPhi",diphPhi);

  // random number generator
  RNG = new TRandom(928); // (argument is seed)
};


void EventTree::GetEvent(Int_t i) {
  if(i%10000==0) printf("[+] %.2f%%\n",100*(float)i/((float)ENT));

  chain->GetEntry(i);

  // theta symmetrization tests
  //theta = fabs( fabs(theta-PI/2.0) - PI/2.0 ); // HERMES symmetrization
  //if(Z[qB] > Z[qA]) theta = PI - theta; // Z symmetrization
  //if( RNG->Rndm() > 0.5 ) theta = PI - theta; // coin-flip symmetrization
  //theta = PI-theta; // full theta flip

  // DIS cuts
  cutQ2 = Q2 > 1.0;
  cutW = W > 2.0;
  cutY = y < 0.8;
  cutDIS = cutQ2 && cutW && cutY;



  // diphoton and pi0/BG cuts
  // -- if dihadron does not include pi0s, just set to true; first we set
  //    the booleans' default values to true
  for(int h=0; h<2; h++) {
    cutDiphKinematics[h] = true;
    cutDiph[h] = true;
  };
  // -- then if there are diphotons in this dihadron we evaluate their cuts
  if(diphCnt>0) {
    for(int h=0; h<diphCnt; h++) {
      
      if(runnum < 4500) { // spring 2018 and before (dnp2018 cuts)
        cutDiphKinematics[h] = diphPhotE[h][0]>0.5 && diphPhotE[h][1]>0.5 &&
                               diphAlpha[h] > 0.05 && diphAlpha[h] < 0.2 &&
                               diphPt[h] > 0.15 &&
                               diphZ[h] > 0.1 && diphZ[h] < 0.6 &&
                               Tools::PhiFiducialCut(diphPhotPhi[h][0]) && 
                               Tools::PhiFiducialCut(diphPhotPhi[h][1]);
      } else { // fall 2018 and after
        cutDiphKinematics[h] = diphPhotE[h][0]>0.5 && diphPhotE[h][1]>0.5;
                               //Tools::PhiFiducialCut(diphPhotPhi[h][0]) && 
                               //Tools::PhiFiducialCut(diphPhotPhi[h][1]);
      };

      // mass cut (depends on whether pi0 signal or BG is desired
      if(useDiphBG) {
        // BG cut
        cutDiph[h] = cutDiphKinematics[h] && diphM[h] < 0.1 || diphM[h] > 0.16;
      } else {
        // pi0 cut
        cutDiph[h] = cutDiphKinematics[h] && diphM[h] >= 0.1 && diphM[h] <= 0.16;
      };
    };
  };



  // dihadron kinematics cuts
  cutDihadronKinematics = 
    Z[qA] > 0.1 && Z[qB] > 0.1 &&
    Zpair < 0.95 &&
    Mmiss > 1.05 &&
    xF > 0 &&
    hadP[qA] > 1.0 && hadP[qB] > 1.0;

  // cutDihadron is the full dihadron cut
  cutDihadron = 
    Tools::PairSame(hadIdx[qA],hadIdx[qB],whichHad[qA],whichHad[qB]) &&
    cutDihadronKinematics && 
    cutDiph[qA] && cutDiph[qB];

  
  // cut for doing cross-checks
  // -- tim's cuts
  cutCrossCheck = 
    Tools::PairSame(hadIdx[qA],hadIdx[qB],kPip,kPim) &&
    /*
    particleCnt[kPip]==1 && particleCnt[kPim]==1 &&
    Tools::EMtoP(hadE[qA],PartMass(kPip)) > 1.0 &&
    Tools::EMtoP(hadE[qB],PartMass(kPim)) > 1.0 &&
    Tools::EMtoP(eleE,PartMass(kE)) > 2.0 &&
    */
    Q2>1 && W>2;
    

  // set preferred PhiR definition
  PhiR = PhiRp;

  PhiHR = Tools::AdjAngle( PhiH - PhiR );

  // calculate zeta
  zeta = ( Z[qA] - Z[qB] ) / Zpair;

};


Bool_t EventTree::Valid() {
  return cutDIS && cutDihadron;
};


void EventTree::PrintEvent() {
  printf("[---] Event Info\n");
  printf("  evnum=%d",evnum);
  printf("  runnum=%d",runnum);
  printf("  pairType=0x%x",pairType);
  printf("\n");
  printf("  helicity=%d",helicity);
  printf("  torus=%.1f",torus);
  printf("\n");
  printf("  triggerBits=%lld",triggerBits);
  printf("\n");
  printf("[---] DIS Kinematics\n");
  printf("  x=%.2f",x);
  printf("  Q2=%.2f",Q2);
  printf("  W=%.2f",W);
  printf("\n");
  printf("  Nu=%.2f",Nu);
  printf("  y=%.2f",y);
  printf("\n");
  printf("[---] Hadron Kinematics: %s\n",PairName(hadIdx[qA],hadIdx[qB]).Data());
  for(int h=0; h<2; h++) {
    printf(" (%s)\n",PartName(dihHadIdx(hadIdx[qA],hadIdx[qB],h)).Data());
    printf("  E=%.2f",hadE[h]);
    printf("  P=%.2f",hadP[h]);
    printf("  Pt=%.2f",hadPt[h]);
    printf("\n");
    printf("  Eta=%.2f",hadEta[h]);
    printf("  Phi=%.2f",hadPhi[h]);
    printf("\n");
  };
  printf("[---] Dihadron Kinematics\n");
  printf("  Mh=%.2f",Mh);
  printf("  Mmiss=%.2f",Mmiss);
  printf("  xF=%.2f",xF);
  printf("  alpha=%.2f\n",alpha);
  printf("  Zpair=%.2f  Z(a)=%.2f  Z(b)=%.2f\n",Zpair,Z[qA],Z[qB]);
  printf("  PhiH=%.2f\n",PhiH);
  printf("[---] PhiR Tests\n");
  printf("  PhiRq=%.2f",PhiRq);
  printf("  PhiRp=%.2f",PhiRp);
  printf("  PhiRp_r=%.2f",PhiRp_r);
  printf("\n");
};


EventTree::~EventTree() {
};

