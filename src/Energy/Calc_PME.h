#ifndef INC_ENERGY_CALC_PME_H
#define INC_ENERGY_CALC_PME_H
#include "PME_Recip.h"
#include "VDW_LongRange_Correction.h"
#include "../ExclusionArray.h"
#include "../PairList.h"
#include "../PairListEngine_Ewald_LJLR.h"
class AtomMask;
class Box;
class EwaldOptions;
class Frame;
class Topology;
namespace Cpptraj {
namespace Energy {
class Calc_PME {
  public:
    Calc_PME();
    /// Init with Box, EwaldOptions and debug level
    int Init(Box const&, EwaldOptions const&, int);
    int Setup(Topology const&, AtomMask const&);
    int CalcNonbondEnergy(Frame const&, AtomMask const&, double&, double&);

  private:
    PairListEngine_Ewald_LJLR<double> NBengine_;
    PME_Recip Recip_;
    PairList pairList_;
    ExclusionArray Excluded_;
    VDW_LongRange_Correction VDW_LR_; ///< For calculating the long range VDW correction
    Timer t_total_;
};
}
}
#endif
