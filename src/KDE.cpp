#include <cmath>
#include <algorithm> // std::min, std::max
#include "KDE.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#ifdef _OPENMP
#  include <omp.h>
#endif

const double KDE::ONE_OVER_ROOT_TWOPI = 1.0 / sqrt( Constants::TWOPI );

double KDE::GaussianKernel(double u) const {
  return ( ONE_OVER_ROOT_TWOPI * exp( -0.5 * u * u ) );
}

/** CONSTRUCTOR - default to Gaussian kernel. */
KDE::KDE() :
# ifdef _OPENMP
  numthreads_(1),
# endif
  ktype_(GAUSSIAN), Kernel_(&KDE::GaussianKernel) {}

#ifdef _OPENMP
KDE::KDE(int n) : numthreads_(n), ktype_(GAUSSIAN), Kernel_(&KDE::GaussianKernel) {}
#endif

int KDE::CalcKDE(DataSet_double& Out, DataSet_1D const& Pdata) const {
  if (Pdata.Size() < 2) {
    mprinterr("Error: Not enough data for KDE.\n");
    return 1;
  }
  // Automatically determine min, max, step, and bin values.
//  std::vector<double> data;
//  data.reserve( Pdata.Size() );
  double N = 0.0;
  double mean = 0.0;
  double M2 = 0.0;
  double min = Pdata.Dval(0);
  double max = min; 
  for (unsigned int i = 0; i != Pdata.Size(); i++) {
    double x = Pdata.Dval(i);
    min = std::min(min, x);
    max = std::max(max, x);
    N++;
    double delta = x - mean;
    mean += delta / N;
    M2 += delta * (x - mean);
//    data.push_back( x );
  }
  M2 /= (N - 1.0);
  double stdev = sqrt(M2);
  double step = 0.0;
  int bins = (int)sqrt((double)Pdata.Size());
/*
  std::sort(data.begin(), data.end());
  double min = data.front();
  double max = data.back();
  unsigned int upperidx, loweridx;
  if ( (data.size() % 2) == 0 ) {
    // Even number of points. Get Q1 as median of lower and Q3 as median of upper.
    unsigned int halfsize = data.size() / 2;
    loweridx = ((halfsize - 1) / 2);
    upperidx = loweridx + halfsize;
  } else {
    // Odd number of points. Include the median in both halves
    unsigned int lsize = (data.size() + 1) / 2;
    loweridx = ((lsize - 1) / 2);
    unsigned int usize = (data.size() - 1) / 2;
    upperidx = loweridx + usize;
  }
  double Q1 = data[loweridx];
  double Q3 = data[upperidx];
  double step = 2 * ((Q3 - Q1) / pow(data.size(), 1/3));
  int bins = 0;
  mprintf("DEBUG: Q1= %g, Q3= %g, step= %g, min= %g, max= %g, mean= %g, stdev= %g\n",
          Q1, Q3, step, min, max, mean, stdev);
  if (max - min < step) {
    // Would only be 1 bin. Probably noisy.
    mprintf("Warning: Data set is very sparse.\n");
    bins = (int)Pdata.Size() / 10;
    step = 0;
  }
*/
  mprintf("DEBUG: mean= %g, stdev= %g\n", mean, stdev);
  HistBin Xdim;
  if (Xdim.CalcBinsOrStep(min, max, step, bins, Pdata.Meta().Legend()))
    return 1;
  Xdim.PrintHistBin();

  // Automatically determine bandwidth
  double bandwidth = 1.06 * stdev * BandwidthFactor(Pdata.Size());
  mprintf("\tBandwidth: %f\n", bandwidth);

  std::vector<double> Increments(Pdata.Size(), 1.0);

  return CalcKDE(Out, Pdata, Increments, Xdim, bandwidth);
}

double KDE::BandwidthFactor(unsigned int N) {
  return pow( (double)N, (-1.0/5.0) );
}

int KDE::CalcKDE(DataSet_double& Out, DataSet_1D const& Pdata,
                 HistBin const& Xdim, double bandwidth) const
{
  std::vector<double> Increments(Pdata.Size(), 1.0);
  return CalcKDE(Out, Pdata, Increments, Xdim, bandwidth);
}

int KDE::CalcKDE(DataSet_double& Out, DataSet_1D const& Pdata,
                 std::vector<double> const& Increments,
                 HistBin const& Xdim, double bandwidth) const
{
  int inSize = (int)Pdata.Size();
  // Allocate output set, set all to zero.
  Out.Zero( Xdim.Bins() );
  Out.SetDim( Dimension::X, Xdim );
  int outSize = (int)Out.Size();

  int frame, bin;
  double increment, val;
  double total = 0.0;
# ifdef _OPENMP
  int original_num_threads;
# pragma omp parallel
  {
#   pragma omp master
    {
      original_num_threads = omp_get_num_threads();
    }
  }
  // Ensure we only execute with the desired number of threads
  if (numthreads_ < original_num_threads)
    omp_set_num_threads( numthreads_ );
# endif
  // Calculate KDE, loop over input data
# ifdef _OPENMP
  int mythread;
  double **P_thread;
# pragma omp parallel private(frame, bin, val, increment, mythread) reduction(+:total)
  {
    mythread = omp_get_thread_num();
    // Prevent race conditions by giving each thread its own histogram
#   pragma omp master
    {
      P_thread = new double*[ numthreads_ ];
      for (int nt = 0; nt < numthreads_; nt++) {
        P_thread[nt] = new double[ outSize ];
        std::fill(P_thread[nt], P_thread[nt] + outSize, 0.0);
      }
    }
#   pragma omp barrier
#   pragma omp for
# endif
    for (frame = 0; frame < inSize; frame++) {
      val = Pdata.Dval(frame);
      increment = Increments[frame];
      total += increment;
      // Apply kernel across histogram
      for (bin = 0; bin < outSize; bin++)
#       ifdef _OPENMP
        P_thread[mythread][bin] +=
#       else
        Out[bin] +=
#       endif
          (increment * (this->*Kernel_)( (Xdim.Coord(bin) - val) / bandwidth ));
    }
# ifdef _OPENMP
  } // END parallel block
  // Combine results from each thread histogram into Out
  for (int i = 0; i < numthreads_; i++) {
    for (int j = 0; j < outSize; j++)
      Out[j] += P_thread[i][j];
    delete[] P_thread[i];
  }
  delete[] P_thread;
  // Restore original number of threads
  if (original_num_threads != numthreads_)
    omp_set_num_threads( original_num_threads );
# endif
  // Normalize
  for (unsigned int j = 0; j < Out.Size(); j++)
    Out[j] /= (total * bandwidth);
  return 0;
}
