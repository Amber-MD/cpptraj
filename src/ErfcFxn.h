#ifndef INC_ERFCFXN_H
#define INC_ERFCFXN_H
#include "SplineFxnTable.h"
#include "Timer.h"
/// Complimentary error function
class ErfcFxn {
  public:
    ErfcFxn() {}
    /// Fill table with ERFC values using given spacing, begin, and end
    int FillErfcTable(double, double, double);
    /// \return Interpolated erfc value
    double ERFC(double xIn) {
#     ifdef _OPENMP
      return table_.Yval( xIn);
#     else
      t_erfc_.Start();
      double erfcval = table_.Yval(xIn);
      t_erfc_.Stop();
      return erfcval;
#     endif
    }
  private:
    /// ERFC function from sander
    static double erfc_func(double);

    SplineFxnTable table_;
    Timer t_erfc_; ///< For recording how much time is spent doing erfc lookup
};
#endif
