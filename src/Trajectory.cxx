#include "Trajectory.h"

ClassImp(Trajectory)

using namespace std;


Trajectory::Trajectory(Int_t particle_index) {
  debug = false;
  Idx = particle_index;
  Vec = new TLorentzVector();
  Vec->SetXYZM(
    0.0,
    0.0,
    0.0,
    PartMass(Idx)
  );
  printf("Trajectory instantiated\n");
};


void Trajectory::SetMomentum(Float_t px, Float_t py, Float_t pz) {
  Vec->SetXYZM(
    px,
    py,
    pz,
    PartMass(Idx)
  );
};


Trajectory::~Trajectory() {
};

