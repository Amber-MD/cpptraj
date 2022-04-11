#ifndef INC_POTENTIALTERM_ANGLE_H
#define INC_POTENTIALTERM_ANGLE_H
#include "PotentialTerm.h"
#include "ParameterTypes.h"
/// Simple Hooke's law angle term
class PotentialTerm_Angle : public PotentialTerm {
  public:
    PotentialTerm_Angle() : PotentialTerm(ANGLE), angParm_(0), Eang_(0) {}

    int SetupTerm(Topology const&, Box const&, CharMask const&, EnergyArray&);
    void CalcForce(Frame&, CharMask const&) const;
  private:
    void addAngles(AngleArray const&, CharMask const&);

    AngleArray activeAngs_;          ///< Array of angles selected by mask during setup
    AngleParmArray const* angParm_;  ///< Pointer to array containing angle parameters
    double* Eang_;                   ///< Pointer to angle term of energy array.
};
#endif
