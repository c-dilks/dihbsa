#include "Modulation.h"

ClassImp(Modulation)

Modulation::Modulation(Int_t tw_, Int_t l_, Int_t m_,
 Int_t level_, Bool_t enablePW_, Int_t polarization_) {

  // twist, l, m, and level (where level is used if there are additional
  // modulations for a specific set of values {tw,l,m})
  tw = tw_;
  l = l_;
  m = m_;
  lev = level_;

  // if true, enables theta dependence (partial waves) in the form of
  // associated legendre polynomials; by default we leave it turned off
  // so that other programs can choose to turn it on
  enablePW = enablePW_;

  // set polarization of structure function
  polarization = polarization_;

  // validation; will likely crash if any of these errors are thrown
  if( !( tw==0 || tw==2 || tw==3 )) {
    fprintf(stderr,"ERROR: Modulation::Modulation -- bad twist (%d)\n",tw);
  };
  if(l<0 || l>LMAX) {
    fprintf(stderr,"ERROR: Modulation::Modulation -- bad L (%d)\n",l);
  };
  if(TMath::Abs(m) > l) {
    fprintf(stderr,"ERROR: Modulation::Modulation -- bad M (L=%d,M=%d)\n",l,m);
  };
  if(polarization<0 || polarization>=nPOL) {
    fprintf(stderr,"ERROR: Modulation::Modulation -- bad polarization setting\n");
  };


  // build a string for the modulation function which is used as a "base";
  // regexps are used to modify this string into a title, a formula, TF3, etc.
  // -- baseStr azimuthal dependence
  mAbs = TMath::Abs(m);
  if(polarization==kLU) {
    switch(tw) {
      case 0:
        aziStr = "1"; // constant modulation
        break;
      case 2:
        if(m==0) aziStr = "0";
        else aziStr = Form("sin(%d*phiH-%d*phiR)",mAbs,mAbs);
        if(m<0) aziStr = "-"+aziStr; // pull minus sign out front
        break;
      case 3:
        aziStr = Form("sin(%d*phiH+%d*phiR)",1-m,m);
        break;
      default: aziStr = "0";
    };
  }
  else if(polarization==kUU) {
    switch(tw) {
      case 0:
        aziStr = "1"; // constant modulation
        break;
      case 2:
        if(lev==0) { // transverse photon
          if(m==0) aziStr = "1";
          else aziStr = Form("cos(%d*phiH-%d*phiR)",mAbs,mAbs);
        }
        else if(lev==1) { // unpolarized photon
          aziStr = Form("cos(%d*phiH+%d*phiR)",2-m,m);
        }
        else aziStr = "0";
        break;
      case 3:
        aziStr = Form("cos(%d*phiH+%d*phiR)",1-m,m);
        break;
      default: aziStr = "0";
    };
  }
  else aziStr = "0";


  // -- baseStr theta dependence
  // - this is from the partial wave expansion in terms of cos(theta); this follows
  //   formulas 19-21 from arXiv:1408.5721
  // - the associated legendre polynomials are from spherical harmonics
  // - note that theta dependence of |l,m> is equal to |l,-m>
  if(enablePW) {
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
          legStr = "sin(2*theta)";
          break;
        case 2:
          legStr = "pow(sin(theta),2)";
          break;
      };
    };
  } else legStr = "1";

  // -- baseStr concatenate azimuthal and theta dependences
  if(enablePW) baseStr = "("+legStr+")*("+aziStr+")";
  else baseStr = aziStr;

  // -- clean up baseStr, to make it more human-readable
  if(aziStr=="0") baseStr="0";
  Tools::GlobalRegexp(baseStr,TRegexp("1\\*"),"");
  Tools::GlobalRegexp(baseStr,TRegexp("0\\*phi.\\+"),"");
  Tools::GlobalRegexp(baseStr,TRegexp("\\+0\\*phi."),"");
  Tools::GlobalRegexp(baseStr,TRegexp("\\+-"),"-");

  // ----> done building baseStr


  // initialize function
  tf3name = this->ModulationName();
  tf3name.ReplaceAll("mod","modFunc");
  function = new TF3(tf3name,this->Formu(),-PI,PI,-PI,PI,0,PI);
};


// evaluate the modulation for specified values of phiH, phiR, theta
Double_t Modulation::Evaluate(Float_t phiH, Float_t phiR, Float_t theta) {
  return function->Eval(phiH,phiR,theta);
};


// build formula string for TF3
TString Modulation::Formu() {
  formuStr = baseStr;
  Tools::GlobalRegexp(formuStr,TRegexp("sin"),"TMath::Sin");
  Tools::GlobalRegexp(formuStr,TRegexp("cos"),"TMath::Cos");
  Tools::GlobalRegexp(formuStr,TRegexp("pow"),"TMath::Power");
  Tools::GlobalRegexp(formuStr,TRegexp("phiH"),"x");
  Tools::GlobalRegexp(formuStr,TRegexp("phiR"),"y");
  Tools::GlobalRegexp(formuStr,TRegexp("theta"),"z");
  return formuStr;
};


// build formula string for RooFit
TString Modulation::FormuRF() {
  formuStr = baseStr;
  Tools::GlobalRegexp(formuStr,TRegexp("sin"),"TMath::Sin");
  Tools::GlobalRegexp(formuStr,TRegexp("cos"),"TMath::Cos");
  Tools::GlobalRegexp(formuStr,TRegexp("pow"),"TMath::Power");
  Tools::GlobalRegexp(formuStr,TRegexp("phiH"),"rfPhiH");
  Tools::GlobalRegexp(formuStr,TRegexp("phiR"),"rfPhiR");
  Tools::GlobalRegexp(formuStr,TRegexp("theta"),"rfTheta");
  return formuStr;
};



TString Modulation::ModulationTitle() {
  TString retstr = baseStr;
  retstr.ReplaceAll("phi","#phi");
  retstr.ReplaceAll("theta","#theta");
  return retstr;
};
TString Modulation::ModulationName() {
  TString retstr = Form("mod_t%d_l%d_m%s%d_lev%d",tw,l,m<0?"N":"",m,lev);
  return retstr;
};

TString Modulation::StateTitle() {

  TString retstr,polStr,lStr;

  if(tw==0) return "const";
  switch(polarization) {
    case kLU:
      polStr = "LU";
      break;
    case kUU: 
      if(tw==2 && lev==0) polStr = "UU,T";
      else polStr = "UU";
      break;
  };

  lStr = enablePW ? Form("%d",l) : "L";
  
  retstr = Form("|%s,%d>^{tw%d}_{%s}",lStr.Data(),m,tw,polStr.Data());
  return retstr;
};


Modulation::~Modulation() {};
