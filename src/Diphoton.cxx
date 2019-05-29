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


  // energy and momentum
  for(int h=0; h<2; h++) {
    Ephot[h] = (photon[h]->Vec).E();
    momPhoton[h] = (photon[h]->Vec).Vect();
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
  Z = fabs(Ephot[0]-Ephot[1]) / E;

  // invariant mass
  Mgg = vecDiphoton.M();

  // opening angle
  Alpha = TMath::ACos(
    momPhoton[0].Dot(momPhoton[1]) / ( momPhoton[0].Mag() * momPhoton[1].Mag() ) 
  );


  // set Trajectory
  Traj->SetVec(vecDiphoton);


  // set booleans
  validDiphoton = Alpha < 0.19;
  //validDiphoton = Alpha<0.19 && E>2 && Z<0.6;

};


void Diphoton::ResetVars() {
  E = -10000;
  Ephot[0] = -10000;
  Ephot[1] = -10000;
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

