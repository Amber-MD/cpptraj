#ifndef INC_CHARMMPARAM_H
#define INC_CHARMMPARAM_H
#include "ParameterSet.h"
#include "FileName.h"
/// Used to read in CHARMM parameters from CHARMM parameter file.
class CharmmParam {
  public:
    CharmmParam() {}
    int ReadParams(ParameterSet&, FileName const&, int);
  private:
};
#endif
