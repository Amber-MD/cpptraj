#ifndef INC_ENERGY_EWALDCALC_DECOMP_LJCPME_H
#define INC_ENERGY_EWALDCALC_DECOMP_LJCPME_H
#ifdef LIBPME
#include "EwaldCalc_Decomp.h"
#include "PME_Recip.h"
#include "VDW_LongRange_Correction.h"
#include "../PairListEngine_Ewald_Decomp_LJCLR.h"
namespace Cpptraj {
namespace Energy {
/** Energy decomposition for nonbonded calculation using PME and LJC w/ VDW
  * long-range correction.
  * Comple with -DCPPTRAJ_DEBUG_ENEDECOMP for more details on the individual
  * contributions.
  */
class EwaldCalc_Decomp_LJCPME : public EwaldCalc_Decomp {
  public:
    EwaldCalc_Decomp_LJCPME();
    int Init(Box const&, EwaldOptions const&, int);
    int Setup(Topology const&, AtomMask const&); // TODO CharMask?
    int CalcNonbondEnergy(Frame const&, AtomMask const&,
                          PairList const&, ExclusionArray const&,
                          double&, double&);
    void Timing(double) const;
  private:
    PairListEngine_Ewald_Decomp_LJCLR<double> NBengine_;
    PME_Recip Recip_;
    VDW_LongRange_Correction VDW_LR_; ///< For calculating the long range VDW correction
};
}
}
#endif /* LIBPME */
#endif
