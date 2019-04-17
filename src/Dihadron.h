#ifndef Dihadron_
#define Dihadron_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <vector>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TMath.h"

// dihbsa
#include "Constants.h"
#include "Trajectory.h"
#include "DIS.h"




class Dihadron : public TObject
{
  public:
    Dihadron();
    ~Dihadron();

    Bool_t debug;

    void SetEvent(
      Trajectory * trajPlus, Trajectory * trajMinus, DIS * disEv);
    Float_t Mh();

    enum plusminus {hP,hM};
    Trajectory * hadron[2];
    DIS * disKin;
    TLorentzVector vecHad[2];
    TLorentzVector vecPh;
    TLorentzVector vecR;

    Float_t phiR;
    Float_t phiH;
    Float_t z[2];
    Float_t zpair;
    
  private:
    int h;

    TVector3 pQ,pL,pPh,pR,pRperp;
    TVector3 pHad[2];
    TVector3 pHadPerp[2];
    TVector3 crossQL,crossQPh,crossQRperp;

    Float_t sgnH,sgnR;
    Float_t numerH,denomH,numerR,denomR;



  ClassDef(Dihadron,1);
};

#endif
