#ifndef INC_POTENTIALFUNCTION_H
#define INC_POTENTIALFUNCTION_H
#include "PotentialTerm.h"
#include "EnergyArray.h"
#include "CharMask.h"
#include <string>
// Forward declares
class Topology;
class Frame;
/// Hold terms for additive potential.
class PotentialFunction {
  public:
    PotentialFunction() : current_(0), deg_of_freedom_(0) {}

    int AddTerm(PotentialTerm::Type);

    int SetupPotential(Topology const&, std::string const&);
 
    int SetupPotential(Topology const&, CharMask const&);

    int CalculateForce(Frame&);

    EnergyArray const& Energy() const { return earray_; }

    Topology* CurrentTop() const { return current_; }

    int DegreesOfFreedom() const { return deg_of_freedom_; }
  private:
    typedef std::vector<PotentialTerm*> Parray;

    int setupPotential(Topology const&);

    Parray terms_;       ///< Array of potential function terms
    EnergyArray earray_; ///< Array of energy terms
    CharMask mask_;      ///< Selected atoms
    Topology* current_;  ///< Pointer to current topology
    int deg_of_freedom_; ///< Current degrees of freedom.
};
#endif
