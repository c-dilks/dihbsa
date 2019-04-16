#include "Trajectory.h"

ClassImp(Trajectory)

using namespace std;


Trajectory::Trajectory() {
  printf("Trajectory instantiated\n");
  debug = false;

  Vec = new TLorentzVector();
  Idx = -10000;
};


Trajectory::Trajectory(Int_t particle_index) {
  Trajectory();
  Idx = particle_index;
  Vec->SetVectM(TVector3(0.0,0.0,0.0),PartMass(Idx));
};


void Trajectory::SetMomentum(Float_t px, Float_t py, Float_t pz) {
  Vec->SetVect(TVector3(px,py,pz));
};


Trajectory::~Trajectory() {
};

