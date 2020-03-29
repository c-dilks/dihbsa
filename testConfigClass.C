// test new Modulation class

R__LOAD_LIBRARY(DihBsa)
#include "Config.h"
#include "Tools.h"

void testConfigClass() {
  Config * c = new Config();
  c->PrintVars();
};
