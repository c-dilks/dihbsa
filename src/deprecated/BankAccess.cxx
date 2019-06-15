// provides general access to HIPO3 bank entries

#include "BankAccess.h"

namespace clas12 {

  void particle::init(const char *bankName, hipo::reader &r) {
    initBranches(bankName,r);
  }

}
