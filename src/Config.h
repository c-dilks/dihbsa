#ifndef Config_
#define Config_

#include "TSystem.h"
#include "TObject.h"

#include <iostream>
#include <fstream>
#include <algorithm>

class Config : public TObject
{
  public:
    Config();
    ~Config() {};

    void PrintVars();

    // experiment (clas, eic)
    TString Experiment;

    // beam energies
    Float_t EbeamEn,PbeamEn;
    
    // kinematic ranges
    Float_t bdQ2[2];
    Float_t bdW[2];
    Float_t bdMh[2];
    Float_t bdMmiss[2];
    Float_t bdEta[2];
    Float_t bdHadP[2];

  ClassDef(Config,1);
};

#endif
