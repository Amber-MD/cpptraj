#ifndef INC_POTENTIALTERM_H
#define INC_POTENTIALTERM_H
#include <string>
// Forward declares
class Topology;
class CharMask;
class Frame;
/// Abstract base class for a term of a potential function.
class PotentialTerm {
  public:
    enum Type { BOND = 0, NTERMS };
    PotentialTerm(Type t) : type_(t) {}

    virtual int SetupTerm(Topology const&, CharMask const&) = 0;
    virtual void CalcForce(Frame&) const = 0;

  private:
    Type type_;
};
#endif
