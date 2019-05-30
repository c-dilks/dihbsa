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


static Int_t PartColor(int p) {
  switch(p) {
    case kE: return kGray+2;
    case kP: return kAzure;
    case kN: return kAzure+10;
    case kPip: return kBlue;
    case kPim: return kRed;
    case kPi0: return kMagenta;
    case kKp: return kGreen+1;
    case kKm: return kGreen-1;
    case kPhoton: return kOrange;
    default: 
      fprintf(stderr,"ERROR: bad PartColor request\n");
      return kBlack;
  };
};

static TString PartColorName(int p) {
  switch(p) {
    case kE: return "grey";
    case kP: return "darkBlue";
    case kN: return "lightBlue";
    case kPip: return "blue";
    case kPim: return "red";
    case kPi0: return "magenta";
    case kKp: return "lightGreen";
    case kKm: return "darkGreen";
    case kPhoton: return "orange";
    default: 
      fprintf(stderr,"ERROR: bad PartColor request\n");
      return "black";
  };
};


// charge sign constants
enum plusminus {hP,hM};

// pair types
enum pairTypeEnum { pairPM, pairP0, pairM0, nPairType };

static TString pairName(int pair) {
  switch(pair) {
    case pairPM: return "pairPM";
    case pairP0: return "pairP0";
    case pairM0: return "pairM0";
    default:
      fprintf(stderr,"ERROR: bad pairName request\n");
      return "unknown";
  };
};

static TString pairTitle(int pair) {
  switch(pair) {
    case pairPM: return "pi+ pi-";
    case pairP0: return "pi+ pi0";
    case pairM0: return "pi0 pi-";
    default:
      fprintf(stderr,"ERROR: bad pairName request\n");
      return "unknown";
  };
};


static Int_t pmIdx(int pair, int h) {
  if(pair == pairPM) {
    if(h == hP)      return kPip;
    else if(h == hM) return kPim;
  } 
  else if(pair == pairP0) {
    if(h == hP)      return kPip; // higher charge
    else if(h == hM) return kPi0; // lower charge
  }
  else if(pair == pairM0) {
    if(h == hP)      return kPi0; // higher charge
    else if(h == hM) return kPim; // lower charge
  };
  fprintf(stderr,"ERROR: bad pmIdx request\n");
  return -1;
};


static TString pmName(int pair, int h) {
  return PartName(pmIdx(pair,h));
};


static TString pmTitle(int pair, int h) {
  return PartTitle(pmIdx(pair,h));
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
