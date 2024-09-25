#ifndef INC_ENERGY_EWALDCALC_DECOMP_PME_H
#define INC_ENERGY_EWALDCALC_DECOMP_PME_H
#include "PME_Recip.h"
#include "VDW_LongRange_Correction.h"
#include "../ExclusionArray.h"
#include "../PairList.h"
#include "../PairListEngine_Ewald_Decomp_LJLR.h"
class AtomMask;
class Box;
class EwaldOptions;
class Frame;
class Topology;
namespace Cpptraj {
namespace Energy {
class EwaldCalc_Decomp_PME {
  public:
    typedef std::vector<double> Darray;

    EwaldCalc_Decomp_PME();
    int Init(Box const&, EwaldOptions const&, int);
    int Setup(Topology const&, AtomMask const&); // TODO CharMask?
    int CalcDecomposedNonbondEnergy(Frame const&, AtomMask const&, double&, double&,
                                    Darray&, Darray&);
  private:
    PairListEngine_Ewald_Decomp_LJLR<double> NBengine_;
    PME_Recip Recip_;
    PairList pairList_;
    ExclusionArray Excluded_;
    VDW_LongRange_Correction VDW_LR_; ///< For calculating the long range VDW correction
    Timer t_total_;
    Timer t_direct_;
};
}
}
#endif
