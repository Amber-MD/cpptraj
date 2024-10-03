#ifndef INC_ENERGY_ECALC_NONBOND_H
#define INC_ENERGY_ECALC_NONBOND_H
#include "../ExclusionArray.h"
#include "../PairList.h"
class EwaldOptions;
class Topology;
namespace Cpptraj {
namespace Energy {
class EwaldCalc;
/// Calculate nonbonded energy for atoms
class Ecalc_Nonbond {
  public:
    enum CalcType { SIMPLE = 0, PME, LJPME, UNSPECIFIED };
    /// CONSTRUCTOR
    Ecalc_Nonbond();
    /// DESTRUCTOR
    ~Ecalc_Nonbond();
    /// Init
    int InitNonbondCalc(CalcType, Box const&, EwaldOptions const&, int);
    /// Setup
    int SetupNonbondCalc(Topology const&, AtomMask const&);
    /// Calculate energy
    int NonbondEnergy(Frame const&, AtomMask const&, double&, double&);
    /// Print timing to stdout
    void PrintTiming(double) const;
  private:
    EwaldCalc* calc_;
    PairList pairList_;
    ExclusionArray Excluded_;
    Timer t_total_;
    CalcType type_;
    bool needs_pairlist_;
};
}
}
#endif
