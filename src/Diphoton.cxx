#include "Diphoton.h"

ClassImp(Diphoton)

using namespace std;


Diphoton::Diphoton() {
  debug = true;

  Traj = new Trajectory(kPi0);
  vecDiphoton = TLorentzVector(0,0,0,0);

  printf("Diphoton instantiated\n");
};


void Diphoton::SetEvent(Trajectory * traj1, Trajectory * traj2) {

  ResetVars();

  photon[0] = traj1;
  photon[1] = traj2;
  vecDiphoton = photon[0]->Vec + photon[1]->Vec;


  // get photon kinematics
  for(int h=0; h<2; h++) {
    momPhoton[h] = (photon[h]->Vec).Vect();
    photE[h] = (photon[h]->Vec).E();
    photPt[h] = (photon[h]->Vec).Pt();
    photEta[h] = (photon[h]->Vec).Eta();
    photPhi[h] = (photon[h]->Vec).Phi();
  };
  E = vecDiphoton.E();

  // transverse momentum
  Pt = vecDiphoton.Pt(); // IN LAB FRAME, wrt BEAM AXIS (maybe change to w.r.t. q?)


  // prevent dividing by 0 issues
  if(E<=0 || Pt<=0) {
    ResetVars();
    return;
  };


  // eta and phi
  Eta = vecDiphoton.Eta();
  Phi = vecDiphoton.Phi();

  // energy sharing
  Z = fabs(photE[0]-photE[1]) / E;

  // invariant mass
  Mgg = vecDiphoton.M();

  // opening angle
  Alpha = TMath::ACos(
    momPhoton[0].Dot(momPhoton[1]) / ( momPhoton[0].Mag() * momPhoton[1].Mag() ) 
  );


  // set Trajectory
  Traj->SetVec(vecDiphoton);


  // set booleans
  //validDiphoton = Alpha < 0.3;
  validDiphoton = true;


};


void Diphoton::ResetVars() {
  for(int h=0; h<2; h++) {
    photE[h] = -10000;
    photPt[h] = -10000;
    photEta[h] = -10000;
    photPhi[h] = -10000;
  };
  E = -10000;
  Z = -10000;
  Pt = -10000;
  Mgg = -10000;
  Alpha = -10000;
  Eta = -10000;
  Phi = -10000;

  validDiphoton = false;
};




Diphoton::~Diphoton() {
};

