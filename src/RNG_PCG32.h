#ifndef INC_RNG_PCG32_H
#define INC_RNG_PCG32_H
#include "RNG.h"
#ifndef _MSC_VER
#include "pcg_random.hpp"
#endif
namespace Cpptraj {

class RNG_PCG32 : public RNG {
  public:
    RNG_PCG32() {}

    //double Generate();
    unsigned int Number();
  private:
    int SetupRng();
#   ifndef _MSC_VER
    pcg32 rng_;
#   endif
};
}
#endif
