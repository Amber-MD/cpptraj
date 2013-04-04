#include <cmath> // sqrt, fabs
#include "DS_Math.h"
#include "ComplexArray.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD, RADDEG

/// \return true if DataSet is cyclic.
static bool IsTorsionArray( DataSet const& ds ) {
  if (ds.ScalarMode() == DataSet::M_TORSION ||
      ds.ScalarMode() == DataSet::M_PUCKER  ||
      ds.ScalarMode() == DataSet::M_ANGLE     )
    return true;
  return false;
}

/// Return true if set is an atomic type (i.e. int, double, float).
static bool GoodCalcType(DataSet const& ds) {
  if (ds.Type() == DataSet::DOUBLE || 
      ds.Type() == DataSet::FLOAT || 
      ds.Type() == DataSet::INT)
    return true;
  mprinterr("Error: DataSet %s is not a valid type for this calc.\n",
            ds.Name().c_str());
  return false;
}

/** Calculate the average over values in this set (and optionally the
  * standard deviation).
  */
double DS_Math::Avg(DataSet& ds, double* stdev) {
  // Check # values
  int numvalues = ds.Size();
  if ( numvalues < 1 ) {
    if (stdev != 0) *stdev = 0.0;
    return 0.0;
  }
  double avg = 0;
  // Check if this set is a good type
  if ( GoodCalcType(ds) ) {
    if (IsTorsionArray(ds)) {
      // Cyclic torsion average
      double sumy = 0.0;
      double sumx = 0.0;
      for ( int i = 0; i < numvalues; ++i ) {
        double theta = ds.Dval( i ) * DEGRAD;
        sumy += sin( theta );
        sumx += cos( theta );
      }
      avg = atan2(sumy, sumx) * RADDEG;
      // Torsion Stdev
      sumy = 0;
      for ( int i = 0; i < numvalues; ++i) {
        double diff = fabs(avg - ds.Dval( i ));
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
        sum += ds.Dval( i );
      avg = sum / (double)numvalues;
      if (stdev==0) return avg;
      // Stdev
      sum = 0;
      for ( int i = 0; i < numvalues; ++i ) {
        double diff = avg - ds.Dval( i );
        diff *= diff;
        sum += diff;
      }
      sum /= (double)numvalues;
      *stdev = sqrt(sum);
    }
  }
  return avg;
}

double DS_Math::Avg(DataSet& ds) {
  return Avg(ds, 0);
}

/** Return the minimum value in the dataset.  */
double DS_Math::Min(DataSet& ds) {
  // Check # values
  if (ds.Size()==0) return 0;
  double min = 0;
  // Check if this set is a good type
  if ( GoodCalcType(ds) ) {
    min = ds.Dval( 0 );
    for (int i = 1; i < ds.Size(); ++i) {
      double val = ds.Dval( i );
      if (val < min) min = val;
    }
  }
  return min;
}

/** Return the maximum value in the dataset.  */
double DS_Math::Max(DataSet& ds) {
  // Check # values
  if ( ds.Size() == 0 ) return 0;
  double max = 0;
  // Check if this set is a good type
  if ( GoodCalcType(ds) ) {
    max = ds.Dval( 0 );
    for (int i = 1; i < ds.Size(); ++i) {
      double val = ds.Dval( i );
      if (val > max) max = val;
    }
  }
  return max;
}

/** Calculate time correlation between two DataSets.
  * \D1 DataSet to calculate correlation for.
  * \D2 DataSet to calculate correlation to.
  * \Ct DataSet to store time correlation fn, must be DOUBLE.
  * \lagmaxIn Max lag to calculate corr. -1 means use size of dataset.
  * \calccovar If true calculate covariance (devation from avg).
  * \return 0 on success, 1 on error.
  */
int DS_Math::CrossCorr( DataSet& D1, DataSet& D2, DataSet& Ct, int lagmaxIn, 
                        bool calccovar, bool usefft )
{
  int lagmax;
  double ct;
  // Check if D1 and D2 are valid types
  if ( !GoodCalcType(D1) ) return 1;
  if ( !GoodCalcType(D2) ) return 1;
  // Check that D1 and D2 have same # data points.
  int Nelements = D1.Size();
  if (Nelements != D2.Size()) {
    mprinterr("Error: CrossCorr: # elements in dataset %s (%i) not equal to\n", 
              D1.Legend().c_str(), Nelements);
    mprinterr("Error:            # elements in dataset %s (%i)\n", 
              D2.Legend().c_str(), D2.Size());
    return 1;
  }
  if (Nelements < 2) {
    mprinterr("Error: CrossCorr: # elements is less than 2 (%i)\n", Nelements);
    return 1;
  }
  // Check return dataset type
  if ( Ct.Type() != DataSet::DOUBLE ) {
    mprinterr("Internal Error: CrossCorr: Ct must be of type DataSet::DOUBLE.\n");
    return 1;
  }
  // Check if lagmaxIn makes sense. Set default lag to be Nelements 
  // if not specified.
  if (lagmaxIn == -1)
    lagmax = Nelements;
  else if (lagmaxIn > Nelements) {
    mprintf("Warning: CrossCorr [%s][%s]: max lag (%i) > Nelements (%i), setting to Nelements.\n",
            D1.Legend().c_str(), D2.Legend().c_str(), lagmaxIn, Nelements);
    lagmax = Nelements;
  } else
    lagmax = lagmaxIn;
  // If calculating covariance calculate averages
  double avg1 = 0;
  double avg2 = 0;
  if ( calccovar ) {
    avg1 = Avg(D1);
    avg2 = Avg(D2);
  }
  // Calculate correlation
  double norm = 1.0;
  if ( usefft ) {
    // Calc using FFT
    CorrF_FFT pubfft1(Nelements);
    ComplexArray data1 = pubfft1.Array();
    data1.PadWithZero(Nelements);
    for (int i = 0; i < Nelements; ++i)
      data1[i*2] = D1.Dval(i) - avg1;
    if (&D2 == &D1)
      pubfft1.AutoCorr(data1);
    else {
      // Populate second dataset if different
      ComplexArray data2 = pubfft1.Array();
      data2.PadWithZero(Nelements);
      for (int i = 0; i < Nelements; ++i)
        data2[i*2] = D2.Dval(i) - avg2;
      pubfft1.CrossCorr(data1, data2);
    }
    // Put real components of data1 in output DataSet
    norm = 1.0 / fabs( data1[0] );
    for (int i = 0; i < lagmax; ++i) {
      ct = data1[i*2] * norm;
      Ct.Add(i, &ct);
    }
  } else {
    // Direct calc
    for (int lag = 0; lag < lagmax; ++lag) {
      ct = 0;
      int jmax = Nelements - lag;
      for (int j = 0; j < jmax; ++j)
        ct += ((D1.Dval(j) - avg1) * (D2.Dval(j+lag) - avg2));
      if (lag == 0) {
        if (ct != 0)
          norm = fabs( ct );
      }
      ct /= norm;
      Ct.Add(lag, &ct);
    }
  }
  return 0;
}

/** Calculate Pearson product-moment correlation between DataSets.
  * \D1 DataSet to caclculate correlation for.
  * \D2 DataSet to caclulate correlation to.
  * \return Pearson product-moment correlation coefficient.
  */
double DS_Math::CorrCoeff( DataSet& D1, DataSet& D2 ) {
  // Check if D1 and D2 are valid types
  if ( !GoodCalcType(D1) ) return 0;
  if ( !GoodCalcType(D2) ) return 0;
  // Check that D1 and D2 have same # data points.
  int Nelements = D1.Size();
  if (Nelements != D2.Size()) {
    mprinterr("Error: Corr: # elements in dataset %s (%i) not equal to\n",
              D1.Legend().c_str(), Nelements);
    mprinterr("Error:       # elements in dataset %s (%i)\n",
              D2.Legend().c_str(), D2.Size());
    return 0;
  }
  // Calculate averages
  double avg1 = Avg(D1);
  double avg2 = Avg(D2);
  // Calculate average deviations. 
  double sumdiff1_2 = 0.0;
  double sumdiff2_2 = 0.0;
  double corr_coeff = 0.0;
  //mprinterr("DATASETS %s and %s\n", c_str(), D2.c_str());
  for (int i = 0; i < Nelements; i++) {
    double diff1 = D1.Dval(i) - avg1;
    double diff2 = D2.Dval(i) - avg2;
    sumdiff1_2 += (diff1 * diff1);
    sumdiff2_2 += (diff2 * diff2);
    corr_coeff += (diff1 * diff2);
  }
  if (sumdiff1_2 == 0.0 || sumdiff2_2 == 0.0) {
    mprintf("Warning: Corr: %s to %s, Normalization is 0\n",
            D1.Legend().c_str(),  D2.Legend().c_str());
    return 0;
  }
  // Correlation coefficient
  corr_coeff /= ( sqrt( sumdiff1_2 ) * sqrt( sumdiff2_2 ) );
  //mprintf("    CORRELATION COEFFICIENT %6s to %6s IS %10.4f\n",
  //        D1_->c_str(), D2_->c_str(), corr_coeff );
  return corr_coeff;
}
