#include "DIS.h"

ClassImp(DIS)

using namespace std;


DIS::DIS() {
  printf("DIS instantiated\n");
  debug = false;
  speedup = true;

  BeamEn = 10.6; // ?? use run db to get beam en

  vecBeam = TLorentzVector(
    0.0,
    0.0,
    BeamEn,
    BeamEn
  );

  vecTarget = TLorentzVector(
    0.0,
    0.0,
    0.0,
    PartMass(kP)
  );
};


DIS::~DIS() {
};

void DIS::SetBeamEn(Float_t newBeamEn) {
  BeamEn = newBeamEn;
  vecBeam.SetPz(BeamEn);
  vecBeam.SetE(BeamEn);
  return;
};


void DIS::SetElectron(Float_t px, Float_t py, Float_t pz) {
  vecElectron.SetXYZM(px,py,pz,PartMass(kE));
  return;
};


// compute DIS kinematics
void DIS::Analyse() {
  ResetVars();
  vecW = vecBeam + vecTarget - vecElectron;
  W = vecW.M();


  // speedup: only compute x,Q2 if W>2
  if(speedup) { if( W<2.0) return; };

  vecQ = vecBeam - vecElectron;
  Q2 = -1 * vecQ.M2();

  Nu = vecBeam.E() - vecElectron.E();
  X = Q2 / ( 2 * PartMass(kP) * Nu );

  vecBreit = vecQ + 2*X*vecTarget;
  BreitBoost = -1 * vecBreit.BoostVector();


  if(debug) {
    Print();
    BreitPrint();
  };

  return;
};


void DIS::Print() {
  printf("[DIS] Lab Frame:\n");
  printf("beam\t");
  vecBeam.Print();
  printf("target\t");
  vecTarget.Print();
  printf("elec\t");
  vecElectron.Print();
  printf("Q\t");
  vecQ.Print();
  printf("W\t");
  vecW.Print();
}

void DIS::BreitPrint() {
  breitBeam = vecBeam;
  breitTarget = vecTarget;
  breitElectron = vecElectron;
  breitW = vecW;
  breitQ = vecQ;

  breitBeam.Boost(BreitBoost);
  breitTarget.Boost(BreitBoost);
  breitElectron.Boost(BreitBoost);
  breitW.Boost(BreitBoost);
  breitQ.Boost(BreitBoost);

  printf("[DIS] Breit Frame:\n");
  //printf("beam\t");
  //breitBeam.Print();
  printf("breit target\t");
  breitTarget.Print();
  //printf("breit elec\t");
  //breitElectron.Print();
  //printf("breit W\t");
  //breitW.Print();
  printf("breit Q\t\t");
  breitQ.Print();
  printf("breit 2xp\t");
  (2*X*breitTarget).Print();

  // check breit frame properties:
  printf("\n");
  // -- 2xp momentum components equal & opposite of Q momentum
  //    (where p is the target momentum in the breit frame)
  printf("Q+2xp:\t");
  (breitQ + 2*X*breitTarget).Print();
  printf("angle(Q,p) = %f\n",breitQ.Vect().Angle(breitTarget.Vect()));
  // -- q is entirely space-like, so E component should be zero
  printf("q0 = %f\n",breitQ.E());

  printf("\n");
};


void DIS::ResetVars() {
  W = -10000;
  Q2 = -10000;
  Nu = -10000;
  X = -10000;
  return;
};

