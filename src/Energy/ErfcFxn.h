#ifndef INC_ERFCFXN_H
#define INC_ERFCFXN_H
#include "../SplineFxnTable.h"
namespace Cpptraj {
namespace Energy {
/// Complimentary error function
class ErfcFxn {
  public:
    ErfcFxn() {}
    /// Fill table with ERFC values using given spacing, begin, and end
    int FillErfcTable(double, double, double);
    /// \return Interpolated erfc value
    double ErfcInterpolated(double xIn) const { return table_.Yval( xIn); }
    /// ERFC function from sander
    static double erfc_func(double);
  private:
    SplineFxnTable table_;
};
}
}
#endif
