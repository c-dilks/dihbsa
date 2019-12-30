#include "Modulation.h"

ClassImp(Modulation)


Modulation::Modulation() {

  // if true, enables theta dependence (partial waves) in the form of
  // associated legendre polynomials
  enableTheta = true;

  // variables which are used to track if a different value of tw,l,m is requested
  twCurr = (Int_t) UNDEF;
  lCurr = (Int_t) UNDEF;
  mCurr = (Int_t) UNDEF;
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

  if( tw!=twCurr || l!=lCurr || m!=mCurr ) {
    twCurr = tw;
    lCurr = l;
    mCurr = m;
    if(funcCurr) delete funcCurr;
    funcCurr = this->BuildTF3(tw,l,m);
  };

  return funcCurr->Eval(phiH,phiR,theta);
  
};


// return a string for the modulation function which is used as a "base";
// a simple regexp can then be used to modify this string into a title, a formula
// for TF1, or anything else you want
TString Modulation::BaseString(Int_t tw, Int_t l, Int_t m) {
  if(!Validate(tw,l,m)) return "unknown";
  mAbs = TMath::Abs(m);

  // azimuthal dependence
  // - this is from the F_LU structure function's azimuthal modulation 
  //   longitudinally polarized electron beam and unpolarized nucleon target
  // - if you don't want to consider the partial wave expansion, technically
  //   you should only consider m==1 as well
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

  // theta dependence
  // - this is from the partial wave expansion in terms of cos(theta)
  // - the associated legendre polynomials are from spherical harmonics
  // - note that theta dependence of |l,m> is equal to |l,-m>
  if(l==0) legStr = "1";
  else if(l==1) {
    switch(mAbs) {
      case 0:
        legStr = "cos(theta)";
        break;
      case 1:
        legStr = "sin(theta)";
        break;
    };
  } else if(l==2) {
    switch(mAbs) {
      case 0:
        legStr = "0.5*(3*pow(cos(theta),2)-1)";
        break;
      case 1:
        legStr = "3*sin(theta)*cos(theta)";
        break;
      case 2:
        legStr = "3*pow(sin(theta),2)";
        break;
    };
  };

  // concatenate
  if(enableTheta) baseStr = "("+legStr+")*("+aziStr+")";
  else baseStr = aziStr;

  // clean up the expression, to make it more human-readable
  if(tw==2 && m==0) baseStr="0";
  Tools::GlobalRegexp(baseStr,TRegexp("1\\*"),"");
  Tools::GlobalRegexp(baseStr,TRegexp("0\\*phi.\\+"),"");
  Tools::GlobalRegexp(baseStr,TRegexp("\\+0\\*phi."),"");
  Tools::GlobalRegexp(baseStr,TRegexp("\\+-"),"-");

  return baseStr;
};


TString Modulation::BuildTF3formu(Int_t tw, Int_t l, Int_t m) {
  if(!Validate(tw,l,m)) return "unknown";
  tf3str = this->BaseString(tw,l,m);

  Tools::GlobalRegexp(tf3str,TRegexp("sin"),"TMath::Sin");
  Tools::GlobalRegexp(tf3str,TRegexp("cos"),"TMath::Cos");
  Tools::GlobalRegexp(tf3str,TRegexp("pow"),"TMath::Power");

  Tools::GlobalRegexp(tf3str,TRegexp("phiH"),"x");
  Tools::GlobalRegexp(tf3str,TRegexp("phiR"),"y");
  Tools::GlobalRegexp(tf3str,TRegexp("theta"),"z");

  return tf3str;
};

TF3 * Modulation::BuildTF3(Int_t tw, Int_t l, Int_t m) {
  tf3name = Form("modFunc_t%d_l%d_m%s%d",tw,l,m<0?"N":"",m);
  return new TF3(tf3name,this->BuildTF3formu(tw,l,m),-PI,PI,-PI,PI,0,PI);
};

Modulation::~Modulation() {};
