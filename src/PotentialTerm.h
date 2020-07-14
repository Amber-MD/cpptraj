#ifndef INC_POTENTIALTERM_H
#define INC_POTENTIALTERM_H
#include <string>
// Forward declares
class Topology;
class CharMask;
class Frame;
class EnergyArray;
class Box;
class MdOpts;
/// Abstract base class for a term of a potential function.
class PotentialTerm {
  public:
    enum Type { BOND = 0, OPENMM, NTERMS };
    PotentialTerm(Type t) : type_(t) {}
    virtual ~PotentialTerm() {} // Virtual since this is inherited.

    virtual int InitTerm(MdOpts const&) { return 0; } // TODO pure virtual?
    virtual int SetupTerm(Topology const&, Box const&, CharMask const&, EnergyArray&) = 0;
    virtual void CalcForce(Frame&, CharMask const&) const = 0;

  private:
    Type type_;
};
#endif
