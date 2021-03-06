#ifndef Modulation_
#define Modulation_

// ROOT
#include "TObject.h"
#include "TString.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TF3.h"

// dihbsa
#include "Constants.h"
#include "Tools.h"


class Modulation : public TObject
{
  public:
    Modulation(Int_t tw_, Int_t l_, Int_t m_,
               Int_t level_ = 0,
               Bool_t enablePW_ = false,
               Int_t polarization_ = kLU);
    ~Modulation();

    Double_t Evaluate(Float_t phiR, Float_t phiH, Float_t theta);
    TString Formu();
    TString FormuRF();

    TString ModulationTitle();
    TString ModulationName();
    TString PolarizationTitle();
    TString StateTitle();
    TString AsymmetryTitle();

    Int_t GetTw() { return tw; };
    Int_t GetL() { return l; };
    Int_t GetM() { return m; };
    TString GetBaseString() { return baseStr; };

    Int_t polarization;
    enum polarization_enum { kLU, kUU, kLL, nPOL };
    static const Int_t LMAX = 2;


  private:
    Double_t azi,leg;
    TString aziStr,legStr,baseStr,formuStr,tf3name;
    Int_t mAbs;
    Int_t tw,l,m,lev,enablePW;
    TF3 * function;


  ClassDef(Modulation,1);
};

#endif
