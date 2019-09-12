// Collection of routines to perform math on 1D datasets.
#include <cmath> // sqrt, fabs
#include "DataSet_1D.h"
#include "Corr.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD, RADDEG
#include "CpptrajFile.h" // Regression output

/** Calculate the average over values in this set (and optionally the
  * standard deviation).
  */
double DataSet_1D::Avg(double* stdev) const {
  // Check # values
  int numvalues = Size();
  if ( numvalues < 1 ) {
    if (stdev!=0) *stdev = 0.0;
    return 0.0;
  }
  double avg = 0;
  if (Meta().IsTorsionArray()) {
    // Cyclic torsion average
    double sumy = 0.0;
    double sumx = 0.0;
    for ( int i = 0; i < numvalues; ++i ) {
      double theta = Dval( i ) * Constants::DEGRAD;
      sumy += sin( theta );
      sumx += cos( theta );
    }
    avg = atan2(sumy, sumx) * Constants::RADDEG;
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
  return avg;
}

/** Return the minimum value in the dataset.  */
double DataSet_1D::Min() const {
  // Check # values
  if (Size()==0) return 0;
  double min = 0;
  min = Dval( 0 );
  for (size_t i = 1; i < Size(); ++i) {
    double val = Dval( i );
    if (val < min) min = val;
  }
  return min;
}

/** Return the maximum value in the dataset.  */
double DataSet_1D::Max() const {
  // Check # values
  if ( Size() == 0 ) return 0;
  double max = 0;
  max = Dval( 0 );
  for (size_t i = 1; i < Size(); ++i) {
    double val = Dval( i );
    if (val > max) max = val;
  }
  return max;
}

static inline double PeriodicDiff(double v1, double v2) {
  double diff = v1 - v2;
  if (diff > 180.0)
    diff = 360.0 - diff;
  else if (diff < -180.0)
    diff = 360.0 + diff;
  return diff;
}

/** Calculate time correlation between two DataSets.
  * \D2 DataSet to calculate correlation to.
  * \Ct DataSet to store time correlation fn, must be DOUBLE.
  * \lagmaxIn Max lag to calculate corr. -1 means use size of dataset.
  * \calccovar If true calculate covariance (devation from avg).
  * \return 0 on success, 1 on error.
  */
int DataSet_1D::CrossCorr( DataSet_1D const& D2, DataSet_1D& Ct,
                           int lagmaxIn, bool calccovar, bool usefft ) const
{
  int lagmax;
  double ct;
  // Check that D1 and D2 have same # data points.
  // TODO: size_t
  int Nelements = (int)Size();
  if (Nelements != (int)D2.Size()) {
    mprinterr("Error: CrossCorr: # elements in dataset %s (%i) not equal to\n",
              legend(), Nelements);
    mprinterr("Error:            # elements in dataset %s (%zu)\n",
              D2.legend(), D2.Size());
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
            legend(), D2.legend(), lagmaxIn, Nelements);
    lagmax = Nelements;
  } else
    lagmax = lagmaxIn;
  // If calculating covariance calculate averages
  double avg1 = 0;
  double avg2 = 0;
  if ( calccovar ) {
    avg1 = this->Avg();
    avg2 = D2.Avg();
  }
  // Calculate correlation
  double norm = 1.0;
  if ( usefft ) {
    // Calc using FFT
    CorrF_FFT pubfft1;
    if (pubfft1.CorrSetup(Nelements)) return 1;
    ComplexArray data1 = pubfft1.Array();
    data1.PadWithZero(Nelements);
    if (Meta().IsTorsionArray()) {
      for (int i = 0; i < Nelements; ++i)
        data1[i*2] = PeriodicDiff(avg1, Dval( i ));
    } else {
      for (int i = 0; i < Nelements; ++i)
        data1[i*2] = Dval(i) - avg1;
    }
    if (&D2 == this)
      pubfft1.AutoCorr(data1);
    else {
      // Populate second dataset if different
      ComplexArray data2 = pubfft1.Array();
      data2.PadWithZero(Nelements);
      if (D2.Meta().IsTorsionArray()) {
        for (int i = 0; i < Nelements; ++i)
          data2[i*2] = PeriodicDiff(avg2, D2.Dval( i ));
      } else {
        for (int i = 0; i < Nelements; ++i)
          data2[i*2] = D2.Dval(i) - avg2;
      }
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
    double diff1, diff2;
    for (int lag = 0; lag < lagmax; ++lag) {
      ct = 0;
      int jmax = Nelements - lag;
      for (int j = 0; j < jmax; ++j) {
        if (Meta().IsTorsionArray())
          diff1 = PeriodicDiff(Dval(j), avg1);
        else
          diff1 = Dval(j) - avg1;
        if (D2.Meta().IsTorsionArray())
          diff2 = PeriodicDiff(D2.Dval(j+lag), avg2);
        else   
          diff2 = D2.Dval(j+lag) - avg2;
        ct += (diff1 * diff2);
      }
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
  * \D2 DataSet to caclulate correlation to.
  * \return Pearson product-moment correlation coefficient.
  */
double DataSet_1D::CorrCoeff( DataSet_1D const& D2 ) const {
  // Check that D1 and D2 have same # data points.
  // TODO: size_t
  int Nelements = (int)Size();
  if (Nelements != (int)D2.Size()) {
    mprinterr("Error: Corr: # elements in dataset %s (%i) not equal to\n",
              legend(), Nelements);
    mprinterr("Error:       # elements in dataset %s (%zu)\n",
              D2.legend(), D2.Size());
    return 0;
  }
  // Calculate averages
  double avg1 = this->Avg();
  double avg2 = D2.Avg();
  // Calculate average deviations. 
  double sumdiff1_2 = 0.0;
  double sumdiff2_2 = 0.0;
  double corr_coeff = 0.0;
  //mprinterr("DATASETS %s and %s\n", c_str(), D2.c_str());
  for (int i = 0; i < Nelements; i++) {
    double diff1 = Dval(i) - avg1;
    double diff2 = D2.Dval(i) - avg2;
    sumdiff1_2 += (diff1 * diff1);
    sumdiff2_2 += (diff2 * diff2);
    corr_coeff += (diff1 * diff2);
  }
  if (sumdiff1_2 == 0.0 || sumdiff2_2 == 0.0) {
    mprintf("Warning: Corr: %s to %s, Normalization is 0\n",
            legend(),  D2.legend());
    return 0;
  }
  // Correlation coefficient
  corr_coeff /= ( sqrt( sumdiff1_2 ) * sqrt( sumdiff2_2 ) );
  //mprintf("    CORRELATION COEFFICIENT %6s to %6s IS %10.4f\n",
  //        D1_->c_str(), D2_->c_str(), corr_coeff );
  return corr_coeff;
}

/** This code (especially the error analysis) was adapted from grace 5.1.22
  * fit.c:linear_regression().
  */
int DataSet_1D::LinearRegression( double& slope, double& intercept,
                                  double& correl, CpptrajFile* outfile ) const
{
  if (Size() < 2) {
    mprinterr("Error: '%s' has less than 2 values, cannot calculate regression.\n",
              legend());
    return 1;
  }
  double mesh_size = (double)Size();
  // Averages
  double xavg = 0.0, yavg = 0.0;
  for (unsigned int i = 0; i < Size(); i++) {
    xavg += Xcrd(i);
    yavg += Dval(i);
  }
  xavg /= mesh_size;
  yavg /= mesh_size;
  // Sums of squares
  double sxx = 0.0, sxy = 0.0, syy = 0.0;
  for (unsigned int i = 0; i < Size(); i++) {
    double xdiff = Xcrd(i) - xavg;
    double ydiff = Dval(i) - yavg;
    sxx += (xdiff * xdiff);
    sxy += (xdiff * ydiff);
    syy += (ydiff * ydiff);
  }
  // Standard deviation, covariance
  double xsd = sqrt( sxx / (mesh_size - 1.0) );
  double ysd = sqrt( syy / (mesh_size - 1.0) );
  if (xsd < Constants::SMALL || ysd < Constants::SMALL) {
    mprinterr("Error: '%s': All values of x or y are the same (SD cannot be zero).\n",
              legend());
    return 1;
  }
  double covariance = sxy / (mesh_size - 1.0);
         correl = covariance / (xsd * ysd);
         slope = sxy / sxx;
         intercept = yavg - slope * xavg;
  if (outfile != 0)
    outfile->Printf("\tData points= %zu\n"
                    "\t<X>= %g\n\t<Y>= %g\n"
                    "\tSDx= %g\n\tSDy= %g\n"
                    "\tCorrelation coefficient= %g\n"
                    "\tSlope= %g\n", Size(),
                    xavg, yavg, xsd, ysd, correl, slope);
  // Case N==2, no need for error analysis.
  if (Size() == 2) {
    slope = (Dval(1) - Dval(0)) / (Xcrd(1) - Xcrd(0));
    intercept = Dval(0) - slope * Xcrd(0);
    if (outfile != 0) outfile->Printf("\tIntercept= %g\n", intercept);
    return 0;
  }
  // Error analysis
  double residualSumSq = syy - slope * sxy;
  double residualMeanSq = residualSumSq / (mesh_size - 2.0);
  //double stdErrRegression = sqrt( residualMeanSq );
  double stdErrIntercept = sqrt( residualMeanSq * (1.0 / mesh_size + xavg * xavg / sxx) );
  double stdErrSlope = sqrt( residualMeanSq / sxx );
  double sumSqRegression = syy - residualSumSq;
  double Fval = sumSqRegression / residualMeanSq;
  //double R2 = sumSqRegression / syy;
  if (outfile != 0) {
    outfile->Printf("\tStandard error of slope= %g\n"
                    "\tt - value for slope= %g\n"
                    "\tIntercept= %g\n"
                    "\tStandard Error of intercept= %g\n"
                    "\tt - value for intercept= %g\n",
                    stdErrSlope, slope / stdErrSlope,
                    intercept, stdErrIntercept, intercept / stdErrIntercept);
    outfile->Printf("\tEquation: Y = %g * X + %g\n", slope, intercept);
    outfile->Printf("\tVariance analysis:\n\t%-10s %5s %14s %14s %14s\n",
            "Source", "d.f", "Sum of squares", "Mean square", "F");
    outfile->Printf("\t%-10s %5d %14.7g %14.7g %14.7g\n", "Regression",
            1, sumSqRegression, sumSqRegression, Fval);
    outfile->Printf("\t%-10s %5zu %14.7g %14.7g\n", "Residual",
            Size() - 2, residualSumSq, residualMeanSq);
    outfile->Printf("\t%-10s %5zu %14.7g\n", "Total",  Size() - 1, syy);
  }
  return 0;
}

// ----- Integration routines --------------------------------------------------
/* Just integration */
double DataSet_1D::Integrate(IntegrationType itype) const {
  double sum = 0.0;
  if (Size() < 2) return 0;
  if (itype == TRAPEZOID) {
    for (unsigned int i = 1; i != Size(); i++) {
      double b_minus_a = Xcrd(i) - Xcrd(i-1);
      sum += (b_minus_a * (Dval(i-1) + Dval(i)) * 0.5);
    }
  }
  return sum;
}

/** Integration with cumulative sum. */
double DataSet_1D::Integrate(IntegrationType itype, std::vector<double>& sumOut) const {
  sumOut.clear();
  double sum = 0.0;
  if (Size() < 2) return 0;
  sumOut.reserve( Size() );
  sumOut.push_back( 0 );
  if (itype == TRAPEZOID) {
    for (unsigned int i = 1; i != Size(); i++) {
      double b_minus_a = Xcrd(i) - Xcrd(i-1);
      sum += (b_minus_a * (Dval(i-1) + Dval(i)) * 0.5);
      sumOut.push_back( sum );
    }
  }
  return sum;
}

