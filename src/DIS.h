#ifndef DIS_
#define DIS_

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

// dihbsa
#include "Constants.h"



class DIS : public TObject
{
  public:
    DIS();
    ~DIS();

    void SetBeamEn(Float_t newBeamEn);
    void SetElectron(Float_t px, Float_t py, Float_t pz);
    Bool_t Analyse();

    Float_t BeamEn;
    Float_t W,Q2,Nu,X;
    
  private:
    TLorentzVector * vecBeam;
    TLorentzVector * vecTarget;
    TLorentzVector * vecElectron;
    TLorentzVector * vecW;
    TLorentzVector * vecQ;


  ClassDef(DIS,1);
};

#endif
