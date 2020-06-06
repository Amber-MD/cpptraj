#ifndef INC_SPLINEFXNTABLE_H
#define INC_SPLINEFXNTABLE_H
#include <vector>
//#inc lude "CpptrajStdio.h" // DEBUG
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
      double Xval = xIn - Xmin_;
      int xidx = ((int)(one_over_Dx_ * Xval));
      double dx = Xval - ((double)xidx * Dx_);
      xidx *= 4;
      // DEBUG
      //if (xidx < 0 || xidx >= (int)table_.size())
      //  mprinterr("Error: index %i out of range for X val %g\n", xidx, xIn);
      //return table_[xidx] + 
      //       dx*(table_[xidx+1] + dx*(table_[xidx+2] + dx*table_[xidx+3]));
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
