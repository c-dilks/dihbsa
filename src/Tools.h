#ifndef TOOLS_H_GUARD
#define TOOLS_H_GUARD

#include "TString.h"
#include "TMath.h"

class Tools {
  public:

    static Float_t AdjAngle(Float_t ang) {
      while(ang>PI) ang-=2*PI;
      while(ang<-PI) ang+=2*PI;
      return ang;
    };

    ClassDef(Tools,1);
};

#endif
