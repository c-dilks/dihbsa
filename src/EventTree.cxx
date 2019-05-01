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
  printf("number of entries: %d\n",ENT);


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
  chain->SetBranchAddress("Z",Z);
  chain->SetBranchAddress("Zpair",&Zpair);
  chain->SetBranchAddress("PhiH",&PhiH);
  chain->SetBranchAddress("PhiR",&PhiR);
  chain->SetBranchAddress("Mmiss",&Mmiss);
  chain->SetBranchAddress("xF",&xF);
  chain->SetBranchAddress("Ph",&Ph);
  chain->SetBranchAddress("Pht",&Pht);

  chain->SetBranchAddress("runnum",&runnum);
  chain->SetBranchAddress("evnum",&evnum);
  chain->SetBranchAddress("helicity",&helicity);
  chain->SetBranchAddress("torus",&torus);
  chain->SetBranchAddress("triggerBits",&triggerBits);

  chain->SetBranchAddress("PhiR_T_byKt",&PhiR_T_byKt);
  chain->SetBranchAddress("PhiR_T_byRej",&PhiR_T_byRej);
  chain->SetBranchAddress("PhiR_Perp",&PhiR_Perp);
  chain->SetBranchAddress("PhiR_byPh",&PhiR_byPh);
  chain->SetBranchAddress("PhiR_byPhad",PhiR_byPhad);
  chain->SetBranchAddress("PhiP1P2",&PhiP1P2);

  chain->SetBranchAddress("b_PhiR_T_byKt",&b_PhiR_T_byKt);
  chain->SetBranchAddress("b_PhiR_T_byRej",&b_PhiR_T_byRej);
  chain->SetBranchAddress("b_PhiR_Perp",&b_PhiR_Perp);
  chain->SetBranchAddress("b_PhiR_byPh",&b_PhiR_byPh);
  chain->SetBranchAddress("b_PhiR_byPhad",b_PhiR_byPhad);
  chain->SetBranchAddress("b_PhiP1P2",&b_PhiP1P2);

};


void EventTree::GetEvent(Int_t i) {
  if(i%10000==0) printf("[+] %.2f%%\n",100*(float)i/((float)ENT));

  chain->GetEntry(i);

  // DIS cuts
  cutQ2 = Q2 > 1.0;
  cutW = W > 2.0;
  cutY = y < 0.8;

  // Dihadron kinematics cuts
  cutDihadron = true;
  cutDihadron = cutDihadron && Z[hP] > 0.1 && Z[hM] > 0.1;
  cutDihadron = cutDihadron && Zpair < 0.95;
  cutDihadron = cutDihadron && Mmiss > 1.05;
  cutDihadron = cutDihadron && xF > 0;
  cutDihadron = cutDihadron && hadP[hP] > 1.0 && hadP[hM] > 1.0;
};


void EventTree::Print() {
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
  printf("\n");
  printf("  Zpair=%.2f  Z[%s]=%.2f  Z[%s]=%.2f",
      Zpair,
      PMsym(hP).Data(),Z[hP],
      PMsym(hM).Data(),Z[hM]);
  printf("\n");
  printf("  PhiH=%.2f",PhiH);
  printf("  PhiR=%.2f",PhiR);
  printf("\n");
  printf("[---] PhiR Tests\n");
  printf("  PhiR_T_byKt=%.2f",PhiR_T_byKt);
  printf("  PhiR_T_byRej=%.2f",PhiR_T_byRej);
  printf("  PhiR_Perp=%.2f",PhiR_Perp);
  printf("  PhiR_byPh=%.2f",PhiR_byPh);
  printf("\n");
  printf("  PhiR_byPhad[%s]=%.2f",PMsym(hP).Data(),PhiR_byPhad[hP]);
  printf("  PhiR_byPhad[%s]=%.2f",PMsym(hM).Data(),PhiR_byPhad[hM]);
  printf("  PhiP1P2=%.2f",PhiP1P2);
  printf("\n");
};


EventTree::~EventTree() {
};

