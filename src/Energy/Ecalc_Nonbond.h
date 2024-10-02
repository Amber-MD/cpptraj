#ifndef INC_ENERGY_ECALC_NONBOND_H
#define INC_ENERGY_ECALC_NONBOND_H
#include "../ExclusionArray.h"
#include "../PairList.h"
class EwaldOptions;
class Topology;
namespace Cpptraj {
namespace Energy {
/// Calculate nonbonded energy for atoms
class Ecalc_Nonbond {
  public:
    enum CalcType { SIMPLE = 0, PME, LJPME, UNSPECIFIED };
    /// CONSTRUCTOR
    Ecalc_Nonbond();
    /// Init
    int InitNonbondCalc(CalcType, Box const&, EwaldOptions const&, int);
    /// Setup
    int SetupNonbondCalc(Topology const&, AtomMask const&);
    /// Calculate energy
    int NonbondEnergy(Frame const&, AtomMask const&, double&, double&);
  private:
    PairList pairList_;
    ExclusionArray Excluded_;
    CalcType type_;
    bool needs_pairlist_;
};
}
}
#endif
