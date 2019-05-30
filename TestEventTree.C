R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"

void TestEventTree(TString dir="outroot") {

  EventTree * ev = new EventTree(TString(dir+"/*.root"),1);

  printf("begin loop through %lld events...\n",ev->ENT);
  for(int i=0; i<ev->ENT; i++) {

    ev->GetEvent(i);
  };
};

