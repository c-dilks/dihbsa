#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD

#include "TString.h"
#include "TMath.h"

// pi
// ---------------------------------------------------
static Double_t PI = TMath::Pi();
static Double_t PIe = TMath::Pi() + 0.3;

// particles constants
// ---------------------------------------------------
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

static Int_t PartCharge(int p) {
  switch(p) {
    case kE: return -1;
    case kP: return 1;
    case kN: return 0;
    case kPip: return 1;
    case kPim: return -1;
    case kPi0: return 0;
    case kKp: return 1;
    case kKm: return -1;
    case kPhoton: return 0;
    default: 
      fprintf(stderr,"ERROR: bad PartCharge request\n");
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


// NEW CODE--------------

enum pair_enum { qH, qL }; // qH = high charge,  qL = low charge

// return the hadron particleIndex within the dihadron pair, where "idx"
// represents either the first or second hadron; by convention, the first
// hadron has higher charge than the second
static Int_t HadIdx(Int_t p1, Int_t p2, Int_t idx) {
  if(idx == qH)      { return PartCharge(idx1) >= partCharge(idx2) ? idx1 : idx2; }
  else if(idx == qL) { return PartCharge(idx1) <  partCharge(idx2) ? idx1 : idx2; }
  else {
    fprintf(stderr,"ERROR: bad HadIdx request\n");
    return -10000;
  };
};

static TString PairName(Int_t p1, Int_t p2) {
  return TString( PartName(HadIdx(p1,p2,qH)) + "_" + PartName(HadIdx(p1,p2,qL)) );
};
static TString PairTitle(Int_t p1, Int_t p2) {
  return TString(
    "(" + PartTitle(HadIdx(p1,p2,qH)) + "," + PartTitle(HadIdx(p1,p2,qL)) + ")"
  );
};

// use these methods to change any pi0 names/titles to BG titles, when looking at BG
static void TransformNameBG(TString & str) {
  str.ReplaceAll(PartName(kPi0),"diphBG");
};

static void TransformTitleBG(TString & str) {
  str.ReplaceAll(PartTitle(kPi0),"#gamma#gamma_{BG}");
};

// return true if pair (a1,a2) is equal to (b1,b2), regardless of order
// -- maybe move this to TOOLS; to be used to see if hadron types a1 and a2
//    equal the one we want to look at, set by b1 and b2
static Bool_t PairSame(Int_t a1, Int_t a2, Int_t b1, Int_t b2) {
  return (a1==b1 && a2==b2) || (a1==b2 && a2==b1);
};




// pair types
// OLD CODE--------------
// ---------------------------------------------------
enum plusminus {hP,hM};
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
// ---------------------------------------------------
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
