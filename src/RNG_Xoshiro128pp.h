#ifndef INC_RNG_XOSHIRO128PP_H
#define INC_RNG_XOSHIRO128PP_H
#include "RNG.h"
namespace Cpptraj {

class RNG_Xoshiro128pp : public RNG {
  public:
    RNG_Xoshiro128pp() {}

    //double Generate();
    unsigned int Number();
  private:
    int SetupRng();
};
}
#endif
