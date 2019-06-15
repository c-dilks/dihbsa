// provides general access to HIPO3 bank entries

#ifndef BankAccess_H
#define BankAccess_H

#include "bank.h"
#include "vectors.h"

namespace clas12 {

  class BankAccess : public hipo::bank {

    public:

      BankAccess() = default;
      BankAccess(const char *bankName, hipo::reader &r) : hipo::bank(bankName,r) {};
      ~BankAccess() = default;

      void init(const char *bankName, hipo::reader &r);

      int GetBankEntryOrder(const char * ent) { return getEntryOrder(ent); };
  }
}

#endif
