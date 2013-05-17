#include <cmath> // sqrt, fabs
#include "DataSet_1D.h"
#include "Corr.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD, RADDEG
/// Collection of routines to perform math on 1D datasets.

/// \return true if DataSet is cyclic.
bool DataSet_1D::IsTorsionArray( DataSet_1D const& ds ) {
  if (ds.ScalarMode() == DataSet::M_TORSION ||
      ds.ScalarMode() == DataSet::M_PUCKER  ||
      ds.ScalarMode() == DataSet::M_ANGLE     )
    return true;
  return false;
}

/// Return true if set is an atomic type (i.e. int, double, float).
bool DataSet_1D::GoodCalcType(DataSet_1D const& ds) {
  if (ds.Type() == DataSet::DOUBLE ||
      ds.Type() == DataSet::FLOAT ||
      ds.Type() == DataSet::INTEGER)
    return true;
  mprinterr("Error: DataSet %s is not a valid type for this calc.\n",
            ds.Name().c_str());
  return false;
}

/** Calculate the average over values in this set (and optionally the
  * standard deviation).
  */
double DataSet_1D::Avg(double* stdev) {
  // Check # values
  int numvalues = Size();
  if ( numvalues < 1 ) {
    if (stdev!=0) *stdev = 0.0;
    return 0.0;
  }
  double avg = 0;
  // Check if this set is a good type
  if ( GoodCalcType(*this) ) {
    if (IsTorsionArray(*this)) {
      // Cyclic torsion average
      double sumy = 0.0;
      double sumx = 0.0;
      for ( int i = 0; i < numvalues; ++i ) {
        double theta = Dval( i ) * DEGRAD;
        sumy += sin( theta );
        sumx += cos( theta );
      }
      avg = atan2(sumy, sumx) * RADDEG;
      if (stdev==0) return avg;
      // Torsion Stdev
      sumy = 0;
      for ( int i = 0; i < numvalues; ++i) {
        double diff = fabs(avg - Dval( i ));
        if (diff > 180.0)
          diff = 360.0 - diff;
        diff *= diff;
        sumy += diff;
      }
      sumy /= (double)numvalues;
      *stdev = sqrt(sumy);
    } else {
      // Non-cyclic, normal average
      double sum = 0;
      for ( int i = 0; i < numvalues; ++i )
        sum += Dval( i );
      avg = sum / (double)numvalues;
      if (stdev==0) return avg;
      // Stdev
      sum = 0;
      for ( int i = 0; i < numvalues; ++i ) {
        double diff = avg - Dval( i );
        diff *= diff;
        sum += diff;
      }
      sum /= (double)numvalues;
      *stdev = sqrt(sum);
    }
  }
  return avg;
}

/** Return the minimum value in the dataset.  */
double DataSet_1D::Min() {
  // Check # values
  if (Size()==0) return 0;
  double min = 0;
  // Check if this set is a good type
  if ( GoodCalcType(*this) ) {
    min = Dval( 0 );
    for (size_t i = 1; i < Size(); ++i) {
      double val = Dval( i );
      if (val < min) min = val;
    }
  }
  return min;
}

/** Return the maximum value in the dataset.  */
double DataSet_1D::Max() {
  // Check # values
  if ( Size() == 0 ) return 0;
  double max = 0;
  // Check if this set is a good type
  if ( GoodCalcType(*this) ) {
    max = Dval( 0 );
    for (size_t i = 1; i < Size(); ++i) {
      double val = Dval( i );
      if (val > max) max = val;
    }
  }
  return max;
}

