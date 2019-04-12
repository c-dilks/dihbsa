#ifndef Observable_
#define Observable_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <vector>
#include "TSystem.h"
#include "TObject.h"
#include "TTree.h"

#include "TFile.h"
#include "Constants.h"



class Observable : public TObject
{
  public:
    Observable();
    ~Observable();
    
  private:


  ClassDef(Observable,1);
};

#endif
