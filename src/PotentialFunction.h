#ifndef INC_POTENTIALFUNCTION_H
#define INC_POTENTIALFUNCTION_H
#include "PotentialTerm.h"
#include "EnergyArray.h"
#include "CharMask.h"
#include <string>
// Forward declares
class Topology;
class Frame;
class Box;
class MdOpts;
/// Hold terms for additive potential.
class PotentialFunction {
  public:
    PotentialFunction() : current_(0), deg_of_freedom_(0) {}
    /// Add term to function with given options
    int AddTerm(PotentialTerm::Type, MdOpts const&);
    /// Add term to function with default options
    int AddTerm(PotentialTerm::Type);
    /// Initialize all terms in the potential with the given options
    int InitPotential(MdOpts const&);
    /// Set up all terms in the potential function using mask expression.
    int SetupPotential(Topology const&, Box const&, std::string const&);
    /// Set up all terms in the potential function using given mask
    int SetupPotential(Topology const&, Box const&, CharMask const&);

    int CalculateForce(Frame&);

    void FnInfo() const;

    EnergyArray const& Energy() const { return earray_; }

    Topology* CurrentTop() const { return current_; }

    int DegreesOfFreedom() const { return deg_of_freedom_; }
  private:
    typedef std::vector<PotentialTerm*> Parray;

    int setupPotential(Topology const&, Box const&);

    Parray terms_;       ///< Array of potential function terms
    EnergyArray earray_; ///< Array of energy terms
    CharMask mask_;      ///< Selected atoms
    Topology* current_;  ///< Pointer to current topology
    int deg_of_freedom_; ///< Current degrees of freedom.
};
#endif
