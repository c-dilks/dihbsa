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

static TString PartName(Int_t p) {
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

static TString PartTitle(Int_t p) {
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

static Int_t PartPID(Int_t p) {
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

static Int_t PIDtoIdx(Int_t pid) {
  for(int i=0; i<nParticles; i++) { if(pid==PartPID(i)) return i; };
  return -10000;
};

static Float_t PartMass(Int_t p) {
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

static Int_t PartCharge(Int_t p) {
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

static Int_t PartColor(Int_t p) {
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

static TString PartColorName(Int_t p) {
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
// represents either the first or second hadron (idx==qH or qL, respectively); 
// -- Convention: 
//    - if charges are different: qH has higher charge than qL
//    - if charges are equal:
//      - if particles are different: qH has higher mass than qL
//      - if particles are same: indistinguishable and order doesn't matter
static Int_t dihHadIdx(Int_t p1, Int_t p2, Int_t idx) {
  if(p1==p2) {
    switch(idx) {
      case qH: return p1;
      case qL: return p2;
    };
  } else {

    if(idx==qH) {
      if( PartCharge(p1) == PartCharge(p2) ) {
        return PartMass(p1) >= PartMass(p2) ? p1 : p2;
      } else {
        return PartCharge(p1) > PartCharge(p2) ? p1 : p2;
      };
    }

    else if(idx==qL) {
      if( PartCharge(p1) == PartCharge(p2) ) {
        return PartMass(p1) < PartMass(p2) ? p1 : p2;
      } else {
        return PartCharge(p1) < PartCharge(p2) ? p1 : p2;
      };
    };

  };
  fprintf(stderr,"ERROR: bad dihHadIdx request\n");
  return -10000;
};


static TString PairName(Int_t p1, Int_t p2) {
  return TString( PartName(dihHadIdx(p1,p2,qH)) + "_" + PartName(dihHadIdx(p1,p2,qL)) );
};
static TString PairTitle(Int_t p1, Int_t p2) {
  return TString(
    "(" + PartTitle(dihHadIdx(p1,p2,qH)) + "," + PartTitle(dihHadIdx(p1,p2,qL)) + ")"
  );
};


// enumerator for particles we consider in dihadron pairs (denoted "observable"); 
// useful for looping over pairs we care about
enum observable_enum {
  sPip,
  sPim,
  sPi0,
  sKp,
  sKm,
  nObservables
};
// observable's index ("OI")
static Int_t OI(Int_t s) {
  switch(s) {
    case sPip: return kPip;
    case sPim: return kPim;
    case sPi0: return kPi0;
    case sKp: return kKp;
    case sKm: return kKm;
    default: 
      fprintf(stderr,"ERROR: bad IdxOfSpecies request\n");
      return -10000;
  };
};
static TString ObsName(Int_t s) { return PartName(OI(s)); };
static TString ObsTitle(Int_t s) { return PartTitle(OI(s)); };
  

// use these methods to change any pi0 names/titles to BG titles, when looking at BG
static void TransformNameBG(TString & str) {
  str.ReplaceAll(PartName(kPi0),"diphBG");
};

static void TransformTitleBG(TString & str) {
  str.ReplaceAll(PartTitle(kPi0),"#gamma#gamma_{BG}");
};




// pair types
// OLD CODE--------------
// ---------------------------------------------------
enum plusminus {hP,hM};
enum pairTypeEnum { pairPM, pairP0, pairM0, nPairType };

static TString pairName(Int_t pair) {
  switch(pair) {
    case pairPM: return "pairPM";
    case pairP0: return "pairP0";
    case pairM0: return "pairM0";
    default:
      fprintf(stderr,"ERROR: bad pairName request\n");
      return "unknown";
  };
};

static TString pairTitle(Int_t pair) {
  switch(pair) {
    case pairPM: return "pi+ pi-";
    case pairP0: return "pi+ pi0";
    case pairM0: return "pi0 pi-";
    default:
      fprintf(stderr,"ERROR: bad pairName request\n");
      return "unknown";
  };
};


static Int_t pmIdx(Int_t pair, Int_t h) {
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


static TString pmName(Int_t pair, Int_t h) {
  return PartName(pmIdx(pair,h));
};


static TString pmTitle(Int_t pair, Int_t h) {
  return PartTitle(pmIdx(pair,h));
};


// spin constants
// ---------------------------------------------------
enum spinEnum { sP, sM, nSpin };

static TString SpinName(Int_t s) {
  switch(s) {
    case sP: return "P";
    case sM: return "M";
    default:
      fprintf(stderr,"ERROR: bad SpinName request\n");
      return "unknown";
  };
};

static TString SpinTitle(Int_t s) {
  switch(s) {
    case sP: return "spin +";
    case sM: return "spin -";
    default:
      fprintf(stderr,"ERROR: bad SpinTitle request\n");
      return "unknown";
  };
};


#endif
