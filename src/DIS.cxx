#include "DIS.h"

ClassImp(DIS)

using namespace std;


DIS::DIS() {
  printf("DIS instantiated\n");
  debug = false;

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


// compute DIS kinematics; returns false if the event
// does not pass any DIS cuts (e.g., W<2)
// -- these cuts have comment "disCut"
Bool_t DIS::Analyse() {
  ResetVars();
  vecW = vecBeam + vecTarget - vecElectron;
  W = vecW.M();


  // disCut: make sure event is above elastic & resonance region
  if( W < 2.0 ) return false;

  vecQ = vecBeam - vecElectron;
  Q2 = -1 * vecQ.M2();

  Nu = vecBeam.E() - vecElectron.E();
  X = Q2 / ( 2 * PartMass(kP) * Nu );


  if(debug) {
    printf("beam\n");
    vecBeam.Print();
    printf("target\n");
    vecTarget.Print();
    printf("electron\n");
    vecElectron.Print();
    printf("W\n");
    vecW.Print();
    printf("--> W = %f\n",W);
  };
  

  return true;
};



void DIS::ResetVars() {
  W = -10000;
  Q2 = -10000;
  Nu = -10000;
  X = -10000;
  return;
};

