#ifndef INC_POTENTIALFUNCTION_H
#define INC_POTENTIALFUNCTION_H
#include "PotentialTerm.h"
#include "EnergyArray.h"
// Forward declares
class Topology;
class CharMask;
class Frame;
/// Hold terms for additive potential.
class PotentialFunction {
  public:
    PotentialFunction() {}

    int AddTerm(PotentialTerm::Type);

    int SetupPotential(Topology const&, CharMask const&);

    int CalculateForce(Frame&, CharMask const&);
  private:
    typedef std::vector<PotentialTerm*> Parray;

    Parray terms_;       ///< Array of potential function terms
    EnergyArray earray_; ///< Array of energy terms
};
#endif
