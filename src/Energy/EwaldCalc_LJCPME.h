#ifndef INC_ENERGY_EWALDCALC_LJCPME_H
#define INC_ENERGY_EWALDCALC_LJCPME_H
#ifdef LIBPME
#include "EwaldCalc.h"
#include "PME_Recip.h"
#include "VDW_LongRange_Correction.h"
#include "../PairListEngine_Ewald_LJCLR.h"
namespace Cpptraj {
namespace Energy {
class EwaldCalc_LJCPME : public EwaldCalc {
  public:
    EwaldCalc_LJCPME();
    /// Init with Box, EwaldOptions and debug level
    int Init(Box const&, EwaldOptions const&, int);
    int Setup(Topology const&, AtomMask const&);
    int CalcNonbondEnergy(Frame const&, AtomMask const&,
                          PairList const&, ExclusionArray const&,
                          double&, double&);

    void Timing(double) const;
  private:
    PairListEngine_Ewald_LJCLR<double> NBengine_;
    PME_Recip Recip_;
    VDW_LongRange_Correction VDW_LR_; ///< For calculating the long range VDW correction
};
}
}
#endif
#endif
