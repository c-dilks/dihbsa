#ifndef TOOLS_H_GUARD
#define TOOLS_H_GUARD

#include "TString.h"
#include "TMath.h"
#include "TH1.h"

class Tools {
  public:

    static Float_t AdjAngle(Float_t ang) {
      while(ang>PI) ang-=2*PI;
      while(ang<-PI) ang+=2*PI;
      return ang;
    };

    static Float_t GetFirstFilledX(TH1 * h) {
      for(int b=1; b<=h->GetNbinsX(); b++) {
        if(h->GetBinContent(b)>0) {
          return h->GetBinCenter(b);
        };
      };
      fprintf(stderr,"Tools::GetFirstFilledX called on empty histogram\n");
      return -10000;
    };

    static Float_t GetLastFilledX(TH1 * h) {
      for(int b=h->GetNbinsX(); b>=1; b--) {
        if(h->GetBinContent(b)>0) {
          return h->GetBinCenter(b);
        };
      };
      fprintf(stderr,"Tools::GetLastFilledX called on empty histogram\n");
      return -10000;
    };

    ClassDef(Tools,1);
};

#endif
