#ifndef Modulation_
#define Modulation_

// ROOT
#include "TObject.h"
#include "TString.h"
#include "TMath.h"
#include "TRegexp.h"

// dihbsa
#include "Constants.h"
#include "Tools.h"


class Modulation : public TObject
{
  public:
    Modulation();
    ~Modulation();

    Bool_t Validate(Int_t tw, Int_t l, Int_t m);
    Double_t Evaluate(Int_t tw, Int_t l, Int_t m, 
                      Float_t phiH, Float_t phiR, Float_t theta);
    TString BaseString(Int_t tw, Int_t l, Int_t m);


    Bool_t enableTheta;
    static const Int_t LMAX = 2;

  private:
    Double_t azi,leg;
    TString aziStr,legStr;
    Int_t mAbs;


  ClassDef(Modulation,1);
};

#endif
