#ifndef INC_RNG_MERSENNETWISTER_H
#define INC_RNG_MERSENNETWISTER_H
#include "RNG.h"
#ifdef C11_SUPPORT
#include <random>
#endif
namespace Cpptraj {
class RNG_MersenneTwister : public RNG {
  public:
    RNG_MersenneTwister();

    //double Generate();
    unsigned int Number();
  private:
    int SetupRng();
#   ifdef C11_SUPPORT
    std::mt19937 gen_; ///< The Mersenne twister engine
#   endif
};
}
#endif
