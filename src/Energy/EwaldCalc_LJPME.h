#ifndef INC_ENERGY_CALC_LJPME_H
#define INC_ENERGY_CALC_LJPME_H
#include "PME_Recip.h"
#include "../ExclusionArray.h"
#include "../PairList.h"
#include "../PairListEngine_Ewald_LJPME.h"
namespace Cpptraj {
namespace Energy {
class Calc_LJPME {
  public:
    Calc_LJPME();
    /// Init with Box, EwaldOptions and debug level
    int Init(Box const&, EwaldOptions const&, int);
    int Setup(Topology const&, AtomMask const&);
    int CalcNonbondEnergy(Frame const&, AtomMask const&, double&, double&);

    void Timing(double) const;
  private:
    PairListEngine_Ewald_LJPME<double> NBengine_;
    PME_Recip Recip_;
    PME_Recip LJrecip_;
    PairList pairList_;
    ExclusionArray Excluded_;
    Timer t_total_;
    Timer t_direct_;
};
}
}
#endif
