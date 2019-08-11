// print out banks of any HIPO file

#include <cstdlib>
#include <iostream>

#include "TString.h"

// Clas12Tool
#include "reader.h"
#include "bank.h"
#include "particle.h"
#include "clas12reader.h"


int main(int argc, char** argv) {

#ifdef HIPO_VERSION
  printf("%s compiled with HIPO_VERSION = %d\n",argv[0],HIPO_VERSION);
#else
  fprintf(stderr,"ERROR: HIPO_VERSION preprocessor macro undefined\n");
  exit(0);
#endif


  // ARGUMENTS
  TString infileN;
  if(argc<=1) {
    printf("USAGE: %s [hipo file]\n",argv[0]);
    exit(0);
  };
  if(argc>1) infileN = TString(argv[1]);


  hipo::dictionary factory;


#if HIPO_VERSION == 3

  hipo::reader reader; // HIPO3
  reader.open(infileN.Data());
  reader.readDictionary(factory);

#elif HIPO_VERSION == 4
  clas12::clas12reader c12reader(infileN.Data()); // HIPO4
  c12reader.getReader().readDictionary(factory);

#endif


  factory.show();
};
