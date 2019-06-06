#ifndef TOOLS_H_GUARD
#define TOOLS_H_GUARD

#include "TString.h"
#include "TMath.h"
#include "TH1.h"
#include "TLorentzVector.h"

class Tools {
  public:

    // shift angle to the range [-PI,+PI]
    static Float_t AdjAngle(Float_t ang) {
      while(ang>PI) ang-=2*PI;
      while(ang<-PI) ang+=2*PI;
      return ang;
    };


    // get first filled bin of a histogram
    static Float_t GetFirstFilledX(TH1 * h) {
      for(int b=1; b<=h->GetNbinsX(); b++) {
        if(h->GetBinContent(b)>0) {
          return h->GetBinCenter(b);
        };
      };
      fprintf(stderr,"Tools::GetFirstFilledX called on empty histogram\n");
      return -10000;
    };

    // get last filled bin of a histogram
    static Float_t GetLastFilledX(TH1 * h) {
      for(int b=h->GetNbinsX(); b>=1; b--) {
        if(h->GetBinContent(b)>0) {
          return h->GetBinCenter(b);
        };
      };
      fprintf(stderr,"Tools::GetLastFilledX called on empty histogram\n");
      return -10000;
    };


    // vector projection:
    // returns vA projected onto vB
    static TVector3 Project(TVector3 vA, TVector3 vB) {

      if(fabs(vB.Dot(vB))<0.0001) {
        //fprintf(stderr,"WARNING: Tools::Project to null vector\n");
        return TVector3(0,0,0);
      };

      return vB * ( vA.Dot(vB) / ( vB.Dot(vB) ) );
    };


    // vector rejection: 
    // returns vC projected onto plane transverse to vD
    static TVector3 Reject(TVector3 vC, TVector3 vD) {

      if(fabs(vD.Dot(vD))<0.0001) {
        //fprintf(stderr,"WARNING: Tools::Reject to null vector\n");
        return TVector3(0,0,0);
      };

      return vC - Project(vC,vD);

    };


    // convert whichPair hex number to hadron indices
    // -- whichPair is a 2-digit hex number, where each digit represents 
    //    the particle index of the hadron
    // -- hex is used since each digit ranges 16 values (0-F), allowing space
    //    for 16 particles
    static void DecodeWhichPair(Int_t w, Int_t & ha, Int_t & hb) {
      ha = w & 0xF;
      hb = w>>4 & 0xF;
    };
    // convert hadron indices to whichPair
    static Int_t EncodeWhichPair(Int_t ha, Int_t hb) {
      return (ha<<4) + hb;
    };


    // compare integer pair (a1,a2) to (b1,b2), returns true if they're the same,
    // regardless of order
    static Bool_t PairSame(Int_t a1, Int_t a2, Int_t b1, Int_t b2) {
      return (a1==b1 && a2==b2) || (a1==b2 && a2==b1);
    };



    // azimuthal fiducial cut
    static Bool_t PhiFiducialCut(Float_t phi_) {
      for(int p=0; p<=3; p++) {
        if( fabs( fabs(phi_) - p*PI/3.0) < 0.35 ) return true;
      };
      return false;
    };



    ClassDef(Tools,1);
};

#endif
