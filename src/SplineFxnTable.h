#ifndef INC_SPLINEFXNTABLE_H
#define INC_SPLINEFXNTABLE_H
#include <vector>
//#incl ude "CpptrajStdio.h" // DEBUG
/// Can be used to approximate a function using cubic splines.
class SplineFxnTable {
  public:
    /// CONSTRUCTOR
    SplineFxnTable();
    /// Generic form of the function to approximate.
    typedef double (*FxnType)(double);
    /// Fill the table using given function and spacing, from given min to max with given scale.
    int FillTable(FxnType, double, double, double, double);
    /// Fill the table using given function, mesh size, min, max
    int FillTable(FxnType, int, double, double);
    /// \return Approximated Y value from given X value.
    double Yval(double xIn) const {
      double xval = xIn - Xmin_;
      long int xidx = ((long int)(one_over_Dx_ * xval));
      //mprintf("SPLINEFXNTBL: idx= %li\n", xidx); // DEBUG
      // Delta from index
      double dx = xval - ((double)xidx * Dx_);
      //double newDx = xIn - Xvals_[xidx];
      //mprintf("DEBUG: xidx= %8li  x=%20.10E  x*dx= %20.10E  xval= %20.10E  dx= %20.10E  newDx= %20.10E\n",
      //        xidx, xIn, (double)xidx*Dx_, Xvals_[xidx], dx, newDx);
      // Index into the table
      xidx *= 4;
      // DEBUG
      //if (xidx < 0 || xidx >= (int)table_.size())
      //  mprinterr("Error: index %li out of range (%zu) for X val %g (Xmin= %g 1/dx= %g)\n", xidx, table_.size(), xIn, Xmin_, one_over_Dx_);
#     ifdef SPLINEFXNTABLE_CHECK_RANGE
      // Protect against out of range values
      // NOTE - If args to Yval are chosen carefully, omitting this check speeds things up.
      if (xidx >= (long int)table_.size())
        xidx = table_.size() - 4;
      else if (xidx < 0)
        xidx = 0;
#     endif
      return table_[xidx] + 
             dx*(table_[xidx+1] + dx*(table_[xidx+2] + dx*table_[xidx+3]));
    }
    /// \return Approximated Y value, use internal X table
    double Yval_xtable(double xIn) const {
      long int xidx = (long int)(one_over_Dx_ * (xIn - Xmin_));
      double dx = xIn - Xvals_[xidx];
      xidx *= 4;
      return table_[xidx] + 
             dx*(table_[xidx+1] + dx*(table_[xidx+2] + dx*table_[xidx+3]));
    }
    /// \return More accurate X value via search of internal X table
    double Yval_accurate(double) const;

    unsigned int Nvals() const { return Xvals_.size(); }

    // DEBUG Access to internal table
    std::vector<double> const& InternalTable() const { return table_; }
  private:
    typedef std::vector<double> Darray;

    double Dx_;          ///< Spacing
    double one_over_Dx_; ///< 1 over spacing
    double Xmin_;        ///< Minimum value for which function can be approximated
    double Xmax_;        ///< Maximum value for which function can be approximated
    Darray table_;       ///< Hold Y followed by spline B C D coefficients
    Darray Xvals_;       ///< Hold X values
};
#endif
