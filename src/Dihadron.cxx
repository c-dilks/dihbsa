#include "Dihadron.h"

ClassImp(Dihadron)

using namespace std;


Dihadron::Dihadron() {
  debug = true;

  for(h=0; h<2; h++) vecHad[h] = TLorentzVector(0,0,0,0);
  vecPh = TLorentzVector(0,0,0,0);
  vecR = TLorentzVector(0,0,0,0);

  printf("Dihadron instantiated\n");
};


void Dihadron::SetHadrons(
  Trajectory * trajPlus, 
  Trajectory * trajMinus
) {
  hadron[hP] = trajPlus;
  hadron[hM] = trajMinus;

  for(h=0; h<2; h++) vecHad[h] = *(hadron[h]->Vec);
  vecPh = vecHad[hP] + vecHad[hM];
  vecR = 0.5 * ( vecHad[hP] - vecHad[hM] );
};


Float_t Dihadron::Mh() {
  return vecPh.M();
};




Dihadron::~Dihadron() {
};

