#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD

#include "TString.h"
#include "TMath.h"

static Double_t PI = TMath::Pi();

// particles constants
enum particle_enum {
  kE,
  kP,
  kN,
  kPip,
  kPim,
  kPi0,
  kKp,
  kKm,
  kPhoton,
  nParticles
};

static TString PartName(int p) {
  switch(p) {
    case kE: return "electron";
    case kP: return "proton";
    case kN: return "neutron";
    case kPip: return "piPlus";
    case kPim: return "piMinus";
    case kPi0: return "pi0";
    case kKp: return "KPlus";
    case kKm: return "KMinus";
    case kPhoton: return "photon";
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
    case kPip: return "#pi^{+}";
    case kPim: return "#pi^{-}";
    case kPi0: return "#pi^{0}";
    case kKp: return "K^{+}";
    case kKm: return "K^{-}";
    case kPhoton: return "#gamma";
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
    case kPip: return 211;
    case kPim: return -211;
    case kPi0: return 111;
    case kKp: return 321;
    case kKm: return -321;
    case kPhoton: return 22;
    default: 
      fprintf(stderr,"ERROR: bad PartPID request\n");
      return -10000;
  };
};

static float PartMass(int p) {
  switch(p) {
    case kE: return 0.000511;
    case kP: return 0.938272;
    case kN: return 0.939565;
    case kPip: return 0.139571;
    case kPim: return 0.139571;
    case kPi0: return 0.134977;
    case kKp: return 0.493677;
    case kKm: return 0.493677;
    case kPhoton: return 0.0;
    default: 
      fprintf(stderr,"ERROR: bad PartMass request\n");
      return -10000;
  };
};


// charge sign constants
enum plusminus {hP,hM};

static TString PMstr(int h) {
  switch(h) {
    case hP: return "p";
    case hM: return "m";
    default: 
      fprintf(stderr,"ERROR: bad PMstr request\n");
      return "unknown";
  };
};

static TString PMsym(int h) {
  switch(h) {
    case hP: return "+";
    case hM: return "-";
    default: 
      fprintf(stderr,"ERROR: bad PMsym request\n");
      return "unknown";
  };
};


// spin constants
enum spinEnum { sP, sM, nSpin };

static TString SpinName(int s) {
  switch(s) {
    case sP: return "P";
    case sM: return "M";
    default:
      fprintf(stderr,"ERROR: bad SpinName request\n");
      return "unknown";
  };
};

static TString SpinTitle(int s) {
  switch(s) {
    case sP: return "spin +";
    case sM: return "spin -";
    default:
      fprintf(stderr,"ERROR: bad SpinTitle request\n");
      return "unknown";
  };
};


#endif
