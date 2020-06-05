#ifndef INC_SPLINEFXNTABLE_H
#define INC_SPLINEFXNTABLE_H
#include <vector>
/// Can be used to approximate a function using cubic splines.
class SplineFxnTable {
  public:
    /// CONSTRUCTOR
    SplineFxnTable();
    /// Generic form of the function to approximate.
    typedef double (*FxnType)(double);
    /// Fill the table using given function and spacing from given min to max.
    int FillTable(FxnType, double, double, double);
    /// \return Approximated Y value from given X value.
    double Yval(double xIn) const {
      int xidx = ((int)(one_over_Dx_ * xIn));
      double dx = xIn - ((double)xidx * Dx_);
      xidx *= 4;
      return table_[xidx] + 
             dx*(table_[xidx+1] + dx*(table_[xidx+2] + dx*table_[xidx+3]));
    }
  private:
    typedef std::vector<double> Darray;

    double Dx_;          ///< Spacing
    double one_over_Dx_; ///< 1 over spacing
    double Xmin_;        ///< Minimum value for which function can be approximated
    double Xmax_;        ///< Maximum value for which function can be approximated
    Darray table_;       ///< Hold Y followed by spline B C D coefficients
};
#endif
