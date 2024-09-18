#ifndef INC_ENERGY_EWALDPARAMS_H
#define INC_ENERGY_EWALDPARAMS_H
#include "ErfcFxn.h"
class Box;
namespace Cpptraj {
namespace Energy {
class EwaldParams {
  public:
    EwaldParams();

    int CheckInput(Box const&, int, double, double,
                   double, double, double,
                   double, double);

  private:
    double FindEwaldCoefficient(double, double);

    double ew_coeff_;                  ///< Ewald coefficient
    double lw_coeff_;                  ///< LJ Ewald coefficient
    double switch_width_; ///< Switching window size for LJ switch if active
    double cutoff_;       ///< Direct space cutoff
    double cut2_;         ///< Direct space cutoff squared.
    double cut2_0_;       ///< Direct space cutoff minus switch width, squared.
    double dsumTol_;      ///< Direct space sum tolerance.
    int debug_;

    ErfcFxn erfc_;                ///< Hold spline interpolation for erfc
};
}
}
#endif
