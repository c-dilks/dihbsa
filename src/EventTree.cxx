#include "EventTree.h"

ClassImp(EventTree)

using namespace std;


EventTree::EventTree(TString filelist) {
  printf("EventTree instantiated\n");

  debug = true;

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

  chain->SetBranchAddress("hadE",hadE);
  chain->SetBranchAddress("hadP",hadP);
  chain->SetBranchAddress("hadPt",hadPt);
  chain->SetBranchAddress("hadEta",hadEta);
  chain->SetBranchAddress("hadPhi",hadPhi);

  chain->SetBranchAddress("Mh",&Mh);
  chain->SetBranchAddress("Mmiss",&Mmiss);
  chain->SetBranchAddress("Z",Z);
  chain->SetBranchAddress("Zpair",&Zpair);
  chain->SetBranchAddress("xF",&xF);
  chain->SetBranchAddress("alpha",&alpha);
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

  chain->SetBranchAddress("b_PhiRq",&b_PhiRq);
  chain->SetBranchAddress("b_PhiRp",&b_PhiRp);
  chain->SetBranchAddress("b_PhiRp_r",&b_PhiRp_r);
  chain->SetBranchAddress("b_PhiRp_g",&b_PhiRp_g);

  
  chain->SetBranchAddress("runnum",&runnum);
  chain->SetBranchAddress("evnum",&evnum);
  chain->SetBranchAddress("helicity",&helicity);
  chain->SetBranchAddress("torus",&torus);
  chain->SetBranchAddress("triggerBits",&triggerBits);


};


void EventTree::GetEvent(Int_t i) {
  if(i%10000==0) printf("[+] %.2f%%\n",100*(float)i/((float)ENT));

  chain->GetEntry(i);

  // DIS cuts
  cutQ2 = Q2 > 1.0;
  cutW = W > 2.0;
  cutY = y < 0.8;

  // Dihadron kinematics cuts
  cutDihadron = 
    Z[hP] > 0.1 && Z[hM] > 0.1 &&
    Zpair < 0.95 &&
    Mmiss > 1.05 &&
    xF > 0 &&
    hadP[hP] > 1.0 && hadP[hM] > 1.0;

  // set preferred PhiR definition
  PhiR = PhiRp;
};


void EventTree::PrintEvent() {
  printf("[---] Event Info\n");
  printf("  evnum=%d",evnum);
  printf("  runnum=%d",runnum);
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
  printf("[---] Hadron Kinematics\n");
  for(int h=0; h<2; h++) {
    printf(" (h%s)\n",PMsym(h).Data());
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
  printf("  alpha=%.2f",alpha);
  printf("\n");
  printf("  Zpair=%.2f  Z[%s]=%.2f  Z[%s]=%.2f",
      Zpair,
      PMsym(hP).Data(),Z[hP],
      PMsym(hM).Data(),Z[hM]);
  printf("\n");
  printf("  PhiH=%.2f",PhiH);
  printf("\n");
  printf("[---] PhiR Tests\n");
  printf("  PhiRq=%.2f",PhiRq);
  printf("  PhiRp=%.2f",PhiRp);
  printf("  PhiRp_r=%.2f",PhiRp_r);
  printf("\n");
};


EventTree::~EventTree() {
};

