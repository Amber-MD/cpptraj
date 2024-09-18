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
    /// \return ERFC value of given distance times the Ewald coefficient
    double ErfcEW(double rIn) const { return erfc_.ErfcInterpolated( ew_coeff_*rIn ); }
    /// \return LJ switch fn value
    double Switch_Fn(double rij2) const {
      double cut2_1 = cut2_;
      if (rij2 <= cut2_0_)
        return 1.0;
      else if (rij2 > cut2_1)
        return 0.0;
      else {
        double xoff_m_x = cut2_1 - rij2;
        double fac = 1.0 / (cut2_1 - cut2_0_);
        return (xoff_m_x*xoff_m_x) * (cut2_1 + 2.0*rij2 - 3.0*cut2_0_) * (fac*fac*fac);
      }
    }

    /// \return LJ PME coefficient
    double LW_Coeff() const { return lw_coeff_; }
    /// \return Direct space cutoff (in Ang squared)
    double Cut2() const { return cut2_; }
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
