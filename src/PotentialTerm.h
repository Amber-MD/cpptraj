#ifndef INC_POTENTIALTERM_H
#define INC_POTENTIALTERM_H
#include <string>
// Forward declares
class Topology;
class CharMask;
class Frame;
class EnergyArray;
/// Abstract base class for a term of a potential function.
class PotentialTerm {
  public:
    enum Type { BOND = 0, SIMPLE_LJ_Q, NTERMS };
    PotentialTerm(Type t) : type_(t) {}

    virtual int SetupTerm(Topology const&, CharMask const&, EnergyArray&) = 0;
    virtual void CalcForce(Frame&, CharMask const&) const = 0;

  private:
    Type type_;
};
#endif
