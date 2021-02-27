#ifndef INC_RNG_STDLIB_H
#define INC_RNG_STDLIB_H
#include "RNG.h"
namespace Cpptraj {
/// RNG from stdlib
class RNG_Stdlib : public RNG {
  public:
    RNG_Stdlib();

    //double Generate();
    unsigned int Number();
  private:
    int SetupRng();
};
}
#endif
