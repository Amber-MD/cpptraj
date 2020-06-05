#include "SplineFxnTable.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "Spline.h"
#include "StringRoutines.h"

/** CONSTRUCTOR */
SplineFxnTable::SplineFxnTable() :
  Dx_(0),
  one_over_Dx_(0),
  Xmin_(0),
  Xmax_(0)
{}

/** Fill the spline function table with values from the given function. */
int SplineFxnTable::FillTable(FxnType fxnIn, double dxIn, double minIn, double maxIn)
{
  Dx_ = dxIn;
  if (Dx_ < Constants::SMALL) {
    mprinterr("Error: Spacing for spline table too small or negative.\n");
    return 1;
  }
  one_over_Dx_ = 1.0 / Dx_;

  Xmax_ = maxIn;
  Xmin_ = minIn;
  double width = Xmax_ - Xmin_;
  if (width < Constants::SMALL) {
    mprinterr("Error: Max %g is not larger than min %g\n", Xmax_, Xmin_);
    return 1;
  }

  // Give the width a 1.5x cushion
  unsigned int TableSize = (unsigned int)(one_over_Dx_ * width * 1.5);

  Darray Xvals, Yvals;
  Xvals.reserve( TableSize );
  Yvals.reserve( TableSize );
  // Save X and Y values so we can calc the spline coefficients
  double xval = 0.0;
  for (unsigned int i = 0; i != TableSize; i++) {
    double yval = fxnIn( xval );
    Xvals.push_back( xval );
    Yvals.push_back( yval );
    xval += Dx_;
  }
  Spline cspline;
  cspline.CubicSpline_Coeff(Xvals, Yvals);
  //Xvals.clear();
  // Store values in Spline table
  table_.clear();
  table_.reserve( TableSize * 4 ); // Y B C D
  for (unsigned int i = 0; i != TableSize; i++) {
    table_.push_back( Yvals[i] );
    table_.push_back( cspline.B_coeff()[i] );
    table_.push_back( cspline.C_coeff()[i] );
    table_.push_back( cspline.D_coeff()[i] );
  }
  // Memory saved Y values plus spline B, C, and D coefficient arrays.
  mprintf("\tMemory used by table and splines: %s\n",
          ByteString(table_.size() * sizeof(double), BYTE_DECIMAL).c_str());
  return 0;
}
