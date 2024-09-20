#ifndef INC_ENERGY_EWALDPARAMS_PME_H
#define INC_ENERGY_EWALDPARAMS_PME_H
#include "EwaldParams.h"
class Frame;
namespace Cpptraj {
namespace Energy {
class EwaldParams_PME : public EwaldParams {
  public:
    EwaldParams_PME();
    int InitEwald(Box const&, EwaldOptions const&, int);
    int SetupEwald(Topology const&, AtomMask const&);

    void FillRecipCoords(Frame const&, AtomMask const&);
  private:
    Darray coordsD_; ///< Hold coords for selected atoms (for recip. calc)
    int nfft_[3]; ///< Number of FFT grid points in each direction
    int order_;   ///< PME B spline order
};
}
}
#endif
