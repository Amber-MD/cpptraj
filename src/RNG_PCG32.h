#ifndef INC_RNG_PCG32_H
#define INC_RNG_PCG32_H
#include "RNG.h"
#include "pcg_random.hpp"
namespace Cpptraj {

class RNG_PCG32 : public RNG {
  public:
    RNG_PCG32() {}

    double Generate();
    unsigned int Number();
  private:
    int SetupRng();

    pcg32 rng_;
};
}
#endif
