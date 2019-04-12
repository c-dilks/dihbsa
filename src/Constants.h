#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD

// pdg particle ids
static const int PID_E = 11; // electron
static const int PID_P = 2212; // proton
static const int PID_N = 2112; // neutron
static const int PID_PIP = 211; // pi+
static const int PID_PIM = -211; // pi-
static const int PID_PI0 = 111; // pi0
static const int PID_KP = 321; // K+
static const int PID_KM = -321; // K-
static const int PID_PHOTON = 22; // photon

// pdg particle masses [GeV/c^2]
static const double MASS_E = 0.000511;
static const double MASS_P = 0.938272;
static const double MASS_N = 0.939565;
static const double MASS_PIP = 0.139571;
static const double MASS_PIM = 0.139571;
static const double MASS_PI0 = 0.134977;
static const double MASS_KP = 0.493677;
static const double MASS_KM = 0.493677;
static const double MASS_PHOTON = 0.0;

enum particle_enum {
  kE,
  kP,
  kN,
  nParticles
};


static int arr[nParticles];
//int * fillPID() {
void fillPID() {
  arr[kE] = 11;
  arr[kP] = 2212;
  arr[kN] = 2112;
  //return arr;
};

//const int PID[nParticles] = fillPID();
const int * PID = arr;

#endif
