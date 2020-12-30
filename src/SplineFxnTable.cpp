#include "SplineFxnTable.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "Spline.h"
#include "StringRoutines.h"
#include <cmath>

/** CONSTRUCTOR */
SplineFxnTable::SplineFxnTable() :
  Dx_(0),
  one_over_Dx_(0),
  Xmin_(0),
  Xmax_(0)
{}

/** Print table info to STDOUT. */
void SplineFxnTable::PrintTableInfo(const char* header) const {
  mprintf("%sTable from %g to %g, width is %g table size %zu\n", header, Xmin_, Xmax_, Xmax_-Xmin_, Xvals_.size());
  mprintf("%sSpline table X range is %g to %g\n", header, Xvals_.front(), Xvals_.back());
}

/** Print memory usage to STDOUT. */
void SplineFxnTable::PrintMemUsage(const char* header) const {
  // Memory saved Y values plus spline B, C, and D coefficient arrays.
  mprintf("%sMemory used by table and splines: %s\n", header,
          ByteString((table_.size()+Xvals_.size()+4) * sizeof(double), BYTE_DECIMAL).c_str());
}

/** Fill the spline function table with values from the given function. */
//int SplineFxnTable::FillTable(FxnType fxnIn, double dxIn, double minIn, double maxIn, double scale)
int SplineFxnTable::FillTable(FxnType fxnIn, double dxIn, double minIn, double maxIn)
{
/*
  if (scale < 1.0) {
    mprinterr("Error: Scaling for spline table width %g cannot be less than 1.0x\n", scale);
    return 1;
  }
*/
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
  // Calculate table size
  int ArraySize = ((int)ceil(one_over_Dx_ * width)) + 1;
  int minIdx = 0;

  Xvals_.clear();
  Darray Yvals;
  Xvals_.reserve( ArraySize );
  Yvals.reserve( ArraySize );
  // Save X and Y values so we can calc the spline coefficients
  //double xval = Xmin_;
  for (int i = minIdx; i != ArraySize; i++)
  {
    double xval = Xmin_ + ((double)i * Dx_);
    double yval = fxnIn( xval );
    Xvals_.push_back( xval );
    Yvals.push_back( yval );
  }
  Xmin_ = Xvals_.front();
  Spline cspline;
  cspline.CubicSpline_Coeff(Xvals_, Yvals);
  //Xvals.clear();
  // Store values in Spline table
  table_.clear();
  table_.reserve( Yvals.size() * 4 ); // Y B C D
  for (unsigned int i = 0; i != Yvals.size(); i++) {
    table_.push_back( Yvals[i] );
    table_.push_back( cspline.B_coeff()[i] );
    table_.push_back( cspline.C_coeff()[i] );
    table_.push_back( cspline.D_coeff()[i] );
  }

  return 0;
}

/** Fill the spline function table with values from the given function. */
int SplineFxnTable::FillTable(FxnType fxnIn, int mesh_size, double minIn, double maxIn)
{
  if (mesh_size < 2) {
    mprinterr("Error: Invalid mesh size: %i\n", mesh_size);
    return 1;
  }
  Xmax_ = maxIn;
  Xmin_ = minIn;
  double width = Xmax_ - Xmin_;
  if (width < Constants::SMALL) {
    mprinterr("Error: Max %g is not larger than min %g\n", Xmax_, Xmin_);
    return 1;
  }


  Dx_ = (maxIn - minIn) / (double)mesh_size; 
  if (Dx_ < Constants::SMALL) {
    mprinterr("Error: Spacing for spline table too small or negative.\n");
    return 1;
  }
  one_over_Dx_ = 1.0 / Dx_;

  // Save X and Y values so we can calc the spline coefficients
  Xvals_.clear();
  Xvals_.reserve( mesh_size );
  Darray Yvals;
  Yvals.reserve( mesh_size );
  double s = (minIn + maxIn)/2;
  double d = (maxIn - minIn)/2;
  for (int i = 0; i < mesh_size; i++) {
    Xvals_.push_back( s + d*((double) (2*i + 1 - mesh_size)/(mesh_size - 1)) );
    double yval = fxnIn( Xvals_.back() );
    Yvals.push_back( yval );
  }

  Spline cspline;
  cspline.CubicSpline_Coeff(Xvals_, Yvals);

  // Store values in Spline table
  table_.clear();
  table_.reserve( mesh_size * 4 ); // Y B C D
  for (int i = 0; i != mesh_size; i++) {
    table_.push_back( Yvals[i] );
    table_.push_back( cspline.B_coeff()[i] );
    table_.push_back( cspline.C_coeff()[i] );
    table_.push_back( cspline.D_coeff()[i] );
  }

  return 0;
}

/** Use a binary search to find the nearest x value; will not suffer from
  * round-off issues but is very slow.
  */
double SplineFxnTable::Yval_accurate(double xIn) const {
  // Find X value
  int ind = 0;
  int low = 0;
  int high = Xvals_.size() - 1;

  // Assumptions for Xvals_: > 1 value, monotonic, ascending order
  if (xIn < Xvals_.front())
    ind = -1;
  else if (xIn > Xvals_.back())
    ind = high;
  else {
    while (low <= high) {
      ind = (low + high) / 2;
      if (xIn < Xvals_[ind]) {
        high = ind - 1;
      } else if (xIn > Xvals_[ind + 1]) {
        low = ind + 1;
      } else {
        //return ind;
        break;
      }
    }
  }

  if (ind < 0) {
    ind = 0;
  } else if ((unsigned int)ind > Xvals_.size() - 2) {
    ind = Xvals_.size() - 1;
  }
  double dx = xIn - Xvals_[ind];
  //mprintf("DEBUG: ind= %8i  x=%20.10E  dx= %20.10E\n", ind, xIn, dx);
  int xidx = ind * 4;
  return table_[xidx] + dx*(table_[xidx+1] + dx*(table_[xidx+2] + dx*table_[xidx+3]));
}
