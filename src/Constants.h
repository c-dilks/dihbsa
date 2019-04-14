#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD

#include "TString.h"

enum particle_enum {
  kE,
  kP,
  kN,
  kPIP,
  kPIM,
  kPI0,
  kKP,
  kKM,
  kPHOTON,
  nParticles
};


static TString PartName(int p) {
  switch(p) {
    case kE: return "electron";
    case kP: return "proton";
    case kN: return "neutron";
    case kPIP: return "piPlus";
    case kPIM: return "piMinus";
    case kPI0: return "pi0";
    case kKP: return "KPlus";
    case kKM: return "KMinus";
    case kPHOTON: return "photon";
    default: 
      fprintf(stderr,"ERROR: bad PartName request\n");
      return "unknown";
  };
};


static TString PartTitle(int p) {
  switch(p) {
    case kE: return "e^{-}";
    case kP: return "p";
    case kN: return "n";
    case kPIP: return "#pi^{+}";
    case kPIM: return "#pi^{-}";
    case kPI0: return "#pi^{0}";
    case kKP: return "K^{+}";
    case kKM: return "K^{-}";
    case kPHOTON: return "#gamma";
    default: 
      fprintf(stderr,"ERROR: bad PartTitle request\n");
      return "unknown";
  };
};


static int PartPID(int p) {
  switch(p) {
    case kE: return 11;
    case kP: return 2212;
    case kN: return 2112;
    case kPIP: return 211;
    case kPIM: return -211;
    case kPI0: return 111;
    case kKP: return 321;
    case kKM: return -321;
    case kPHOTON: return 22;
    default: 
      fprintf(stderr,"ERROR: bad PartPID request\n");
      return -10000;
  };
};


static int PartMass(int p) {
  switch(p) {
    case kE: return 0.000511;
    case kP: return 0.938272;
    case kN: return 0.939565;
    case kPIP: return 0.139571;
    case kPIM: return 0.139571;
    case kPI0: return 0.134977;
    case kKP: return 0.493677;
    case kKM: return 0.493677;
    case kPHOTON: return 0.0;
    default: 
      fprintf(stderr,"ERROR: bad PartMass request\n");
      return -10000;
  };
};


#endif
