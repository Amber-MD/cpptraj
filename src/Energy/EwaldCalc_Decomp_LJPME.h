#ifndef INC_ENERGY_EWALDCALC_DECOMP_LJPME_H
#define INC_ENERGY_EWALDCALC_DECOMP_LJPME_H
#ifdef LIBPME
#include "EwaldCalc_Decomp.h"
#include "PME_Recip.h"
#include "../PairListEngine_Ewald_Decomp_LJPME.h"
namespace Cpptraj {
namespace Energy {
class EwaldCalc_Decomp_LJPME : public EwaldCalc_Decomp {
  public:
    EwaldCalc_Decomp_LJPME();
    /// Init with Box, EwaldOptions and debug level
    int Init(Box const&, EwaldOptions const&, int);
    int Setup(Topology const&, AtomMask const&);
    int CalcNonbondEnergy(Frame const&, AtomMask const&,
                          PairList const&, ExclusionArray const&,
                          double&, double&);

    void Timing(double) const;
  private:
    PairListEngine_Ewald_Decomp_LJPME<double> NBengine_;
    PME_Recip Recip_;
    PME_Recip LJrecip_;
};
}
}
#endif
#endif
