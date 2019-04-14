#include "DIS.h"

ClassImp(DIS)

namespace
{
  const Int_t N_CLASS = 5;

  const Float_t Zd = 720;  // distance from IP [cm]
  const Float_t Sd = 3.8; // small cell dimension [cm]
  const Float_t Ld = 5.8; // large cell dimension [cm]

  const Float_t pi0_mass = 0.135; // pi0 mass [GeV]
  const Float_t etm_mass = 0.548; // eta mass [GeV]
  //const Float_t jps_mass = 3.097; // j/psi mass [GeV]
};

DIS::DIS() {
  BeamEn = 10.6; // ?? use run db to get beam en

  vecBeam = new TLorentzVector(
    0.0,
    0.0,
    BeamEn,
    BeamEn
  ); // ?? is initializing as a 4-momentum like this okay?

  vecTarget = new TLorentzVector(
    0.0,
    0.0,
    0.0,
    PartMass(kP)
  );

  vecElectron = new TLorentzVector();
  vecW = new TLorentzVector();
  vecQ = new TLorentzVector();
};


void DIS::SetBeamEn(Float_t newBeamEn) {
  BeamEn = newBeamEn;
  vecBeam->SetPz(BeamEn);
  vecBeam->SetE(BeamEn);
  return;
};


void DIS::SetElectron(Float_t px, Float_t py, Float_t pz) {
  vecElectron->SetVectM(
    TVector3(px,py,pz),PartMass(kE)
  ); // ?? check this definition is okay!
  return;
};


// compute DIS kinematics; returns false if the event
// does not pass any DIS cuts (e.g., W<2)
// -- these cuts have comment "DIScut"
Bool_t DIS::Analyse() {
  ResetVars();
  vecW = vecBeam + vecTarget - vecElectron;
  W = vecW->M();
  

  // DIScut: make sure event is above elastic & resonance region
  if( W < 2.0 ) { 
    ResetVars();
    return false;
  };

  vecQ = vecBeam - vecElectron;
  Q2 = -1 * vecQ->Mag2(); // ?? check this!

  Nu = vecBeam->E() - vecElectron->E();
  X = Q2 / ( 2 * PartMass(kP) * Nu );

  return true;
};



void DIS::ResetVars() {
  W = -10000;
  Q2 = -10000;
  Nu = -10000;
  X = -10000;
  return;
};


DIS::~DIS() {
};
