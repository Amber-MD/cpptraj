#ifndef INC_ENERGY_EWALDPARAMS_LJPME_H
#define INC_ENERGY_EWALDPARAMS_LJPME_H
#include "EwaldParams_PME.h"
namespace Cpptraj {
namespace Energy {
class EwaldParams_LJPME : public EwaldParams_PME {
  public:
    EwaldParams_LJPME();

    int InitEwald(Box const&, EwaldOptions const&, int);
    int SetupEwald(Topology const&, AtomMask const&);

    /// \return LJ Ewald coefficient
    double LJ_EwaldCoeff() const { return lw_coeff_; }

  private:
   double lw_coeff_;                  ///< LJ Ewald coefficient     
};
}
}
#endif
