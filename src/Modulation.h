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
    Modulation();
    ~Modulation();

    Bool_t Validate(Int_t tw, Int_t l, Int_t m);
    Double_t Evaluate(Int_t tw, Int_t l, Int_t m, 
                      Float_t phiH, Float_t phiR, Float_t theta);
    TString BaseString(Int_t tw, Int_t l, Int_t m);
    TString BuildFormu(Int_t tw, Int_t l, Int_t m);
    TString BuildFormuRF(Int_t tw, Int_t l, Int_t m);

    TString ModulationTitle(Int_t tw, Int_t l, Int_t m);
    TString ModulationName(Int_t tw, Int_t l, Int_t m);
    TString StateTitle(Int_t tw, Int_t l, Int_t m);

    Bool_t enablePW;
    static const Int_t LMAX = 2;


  private:
    Double_t azi,leg;
    TString aziStr,legStr,baseStr,formuStr,tf3name;
    Int_t mAbs;
    Int_t twCurr,lCurr,mCurr;
    Int_t twOA,lOA,mOA;
    TF3 * funcCurr;


  ClassDef(Modulation,1);
};

#endif
