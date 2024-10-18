#ifndef INC_ENERGY_EWALDCALC_LJPME_H
#define INC_ENERGY_EWALDCALC_LJPME_H
#ifdef LIBPME
#include "EwaldCalc.h"
#include "PME_Recip.h"
#include "../PairListEngine_Ewald_LJPME.h"
namespace Cpptraj {
namespace Energy {
class EwaldCalc_LJPME : public EwaldCalc {
  public:
    EwaldCalc_LJPME();
    /// Init with Box, EwaldOptions and debug level
    int Init(Box const&, EwaldOptions const&, int);
    int Setup(Topology const&, AtomMask const&);
    int CalcNonbondEnergy(Frame const&, AtomMask const&,
                          PairList const&, ExclusionArray const&,
                          double&, double&);

    void Timing(double) const;
  private:
    PairListEngine_Ewald_LJPME<double> NBengine_;
    PME_Recip Recip_;
    PME_Recip LJrecip_;
};
}
}
#endif
#endif
