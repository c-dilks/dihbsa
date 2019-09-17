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

  // ARGUMENTS
  TString infileN;
  if(argc<=1) {
    printf("USAGE: %s [hipo file]\n",argv[0]);
    exit(0);
  };
  if(argc>1) infileN = TString(argv[1]);


  hipo::dictionary factory;


  clas12::clas12reader c12reader(infileN.Data()); 
  c12reader.getReader().readDictionary(factory);
  factory.show();

  // example of how to read a bank and its rows, within an event loop
  /*
  hipo::event readerEvent;
  hipo::bank mcLund(factory.getSchema("MC::Lund"));
  int nrows,pid;
  float momentum[3];
  while(c12reader.next()==true) {
    c12reader.getReader().read(readerEvent);
    readerEvent.getStructure(mcLund);
    nrows = mcLund.getRows();
    for(int r=0; r<nrows; r++) {
      pid = mcLund.getInt("pid",r);
      momentum[0] = mcLund.getFloat("px",r);
      momentum[1] = mcLund.getFloat("py",r);
      momentum[2] = mcLund.getFloat("pz",r);
    };
  };
  */

};
