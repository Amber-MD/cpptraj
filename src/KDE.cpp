#include "KDE.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include <cmath>

static const double ONE_OVER_ROOT_TWOPI = 1.0 / sqrt( Constants::TWOPI );

static double GaussianKernel(double u) {
  return ( ONE_OVER_ROOT_TWOPI * exp( -0.5 * u * u ) );
}


int KDE::GaussianKDE(DataSet_double& Out, DataSet_1D const& Pdata,
                     std::vector<double> const& Increments,
                     HistBin const& Xdim, double bandwidth)
{
  return KDE::CalcKDE(GaussianKernel, Out, Pdata, Increments, Xdim, bandwidth);
}

int KDE::CalcKDE(FxnPtr Kernel, DataSet_double& Out, DataSet_1D const& Pdata,
                 std::vector<double> const& Increments,
                 HistBin const& Xdim, double bandwidth)
{
  int inSize = (int)Pdata.Size();
  // Allocate output set
  Out.Resize( Xdim.Bins() );
  int outSize = (int)Out.Size();

  int frame, bin;
  double increment;
  double total = 0.0;
# ifdef _OPENMP
  int numthreads;
# pragma omp parallel
  {
#   pragma omp master
    {
      numthreads = omp_get_num_threads();
      mprintf("\tParallelizing calculation with %i threads\n", numthreads);
    }
  }
# endif
    double val;
    // Calculate KDE, loop over input data
#   ifdef _OPENMP
    int mythread;
    double **P_thread;
#   pragma omp parallel private(frame, bin, val, increment, mythread) reduction(+:total)
    {
      mythread = omp_get_thread_num();
      // Prevent race conditions by giving each thread its own histogram
#     pragma omp master
      {
        P_thread = new double*[ numthreads ];
        for (int nt = 0; nt < numthreads; nt++) {
          P_thread[nt] = new double[ outSize ];
          std::fill(P_thread[nt], P_thread[nt] + outSize, 0.0);
        }
      }
#     pragma omp barrier
#     pragma omp for
#   endif
      for (frame = 0; frame < inSize; frame++) {
        val = Pdata.Dval(frame);
        increment = Increments[frame];
        total += increment;
        // Apply kernel across histogram
        for (bin = 0; bin < outSize; bin++)
#         ifdef _OPENMP
          P_thread[mythread][bin] +=
#         else
          Out[bin] +=
#         endif
            (increment * (*Kernel)( (Xdim.Coord(bin) - val) / bandwidth ));
      }
#   ifdef _OPENMP
    } // END parallel block
    // Combine results from each thread histogram into Out
    for (int i = 0; i < numthreads; i++) {
      for (int j = 0; j < outSize; j++)
        Out[j] += P_thread[i][j];
      delete[] P_thread[i];
    }
    delete[] P_thread;
#   endif
  return 0;
}
