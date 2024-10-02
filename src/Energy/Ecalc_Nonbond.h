#ifndef INC_ENERGY_ECALC_NONBOND_H
#define INC_ENERGY_ECALC_NONBOND_H
#include "../ExclusionArray.h"
#include "../PairList.h"
namespace Cpptraj {
namespace Energy {
/// Calculate nonbonded energy for atoms
class Ecalc_Nonbond {
  public:
    enum CalcType { SIMPLE = 0, PME, LJPME };
    /// CONSTRUCTOR
    Ecalc_Nonbond();
    /// Init
    InitNonbondCalc(CalcType, Box const&, EwaldOptions const&, int);
    /// Setup
    SetupNonbondCalc(Topology const&, AtomMask const&);
    /// Calculate energy
    int NonbondEnergy(Frame const&, AtomMask const&, double&, double&);
  private:
    PairList pairList_;
    ExclusionArray Excluded_;
};
}
}
#endif
