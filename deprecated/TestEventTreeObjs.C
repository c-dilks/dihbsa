R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Dihadron.h"
#include "DIS.h"

void TestEventTreeObjs(TString dir="outroot.fall18.some") {

  EventTree * ev = new EventTree(TString(dir+"/*.root"),0x34);

  printf("begin loop through %lld events...\n",ev->ENT);
  for(int i=0; i<ev->ENT; i++) {

    ev->GetEvent(i);
    if(ev->Valid()) {
      //printf("%f %f\n",ev->PhiRp,ev->GetDihadronObj()->PhiRp);
      printf("%f %f\n",ev->W,ev->GetDISObj()->W);
    };

  };
};

