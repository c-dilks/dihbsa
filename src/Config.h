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

    TString Experiment;
    Float_t EbeamEn,PbeamEn;

  ClassDef(Config,1);
};

#endif
