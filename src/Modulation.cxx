#include "Modulation.h"

ClassImp(Modulation)


Modulation::Modulation() {

  // if true, enables theta dependence (partial waves) in the form of
  // associated legendre polynomials
  enableTheta = true;
};


// check for valid values of twist, l, and m
Bool_t Modulation::Validate(Int_t tw, Int_t l, Int_t m) {
  if(tw!=2 && tw!=3) {
    fprintf(stderr,"ERROR: Modulation::Validate -- bad twist (%d)\n",tw);
    return false;
  };
  if(l<0 || l>LMAX) {
    fprintf(stderr,"ERROR: Modulation::Validate -- bad L (%d)\n",l);
    return false;
  };
  if(TMath::Abs(m) > l) {
    fprintf(stderr,"ERROR: Modulation::Validate -- bad M (L=%d,M=%d)\n",l,m);
    return false;
  };
  return true;
};

// evaluate the modulation for specified values of phiH, phiR, theta
Double_t Modulation::Evaluate(Int_t tw, Int_t l, Int_t m, 
                              Float_t phiH, Float_t phiR, Float_t theta) {

  if(!Validate(tw,l,m)) return UNDEF;

  // evaluate azimuthal dependence (uses tw and m)
  // - this is from the F_LU structure function's azimuthal modulation 
  //   longitudinally polarized electron beam and unpolarized nucleon target
  // - if you don't want to consider the partial wave expansion, technically
  //   you should only consider m==1
  switch(tw) {
    case 2:
      azi = TMath::Sin( m * (phiH-phiR) );
      break;
    case 3:
      azi = TMath::Sin( (1-m)*phiH + m*phiR );
      break;
  };

  // evaluate theta dependence (uses l and |m|)
  // - this is from the partial wave expansion in terms of cos(theta)
  // - the associated legendre polynomials are from spherical harmonics
  // - note that theta dependence of |l,m> is equal to |l,-m>
  if(enableTheta) {
    mAbs = TMath::Abs(m);
    if(l==0) leg = 1;
    else if(l==1) {
      switch(mAbs) {
        case 0:
          leg = TMath::Cos(theta);
          break;
        case 1:
          leg = TMath::Sin(theta);
          break;
      };
    } else if(l==2) {
      switch(mAbs) {
        case 0:
          leg = 0.5*( 3*TMath::Power(TMath::Cos(theta),2) - 1 );
          break;
        case 1:
          leg = 3 * TMath::Sin(theta) * TMath::Cos(theta);
          break;
        case 2:
          leg = 3 * TMath::Power(TMath::Sin(theta),2);
          break;
      };
    };
  } else {
    leg=1; // if !enableTheta
  };

  // return the evaluated product of the two functions
  return azi * leg;
};


// return a string for the modulation function which is used as a "base";
// a simple regexp can then be used to modify this string into a title, a formula
// for TF1, or anything else you want
TString Modulation::BaseString(Int_t tw, Int_t l, Int_t m) {
  if(!Validate(tw,l,m)) return "unknown";
  mAbs = TMath::Abs(m);

  // regexps:
  // 1* -> 
  // 0*phi.+ ->
  // +0*phi. ->
  //

  // azimuthal dependence
  switch(tw) {
    case 2:
      if(m==0) aziStr = "0";
      else aziStr = Form("sin(%d*phiH-%d*phiR)",mAbs,mAbs);
      if(m<0) aziStr = "-"+aziStr; // pull minus sign out front
      break;
    case 3:
      aziStr = Form("sin(%d*phiH+%d*phiR)",1-m,m);
      break;
  };

  // strip out factors of 1 and 0
  aziStr(TRegexp("1*")) = "";
  aziStr(TRegexp("0*phi.+")) = "";
  aziStr(TRegexp("+0*phi.")) = "";

  // use this method + regexps to generate TF1s
  // then instead of writing a separate Evaluate method, just use TF1::Eval
};




Modulation::~Modulation() {};
