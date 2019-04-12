#include "Observable.h"

ClassImp(Observable)

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

Observable::Observable() {
  printf("PID_E = %d\n",PID_E);
};

Observable::~Observable() {
};
