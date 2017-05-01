#include <cmath> //sqrt
#include "Bootstrap.h"
#include "CpptrajStdio.h"

// Bootstrap::Init()
int Bootstrap::Init(DataSet_1D* dataIn, int sIn, int nresampleIn, int seedIn, int debugIn)
{
  if (dataIn == 0) {
    mprinterr("Internal Error: Null set passed to Bootstrap::Init\n");
    return 1;
  }
  data_ = dataIn;

  if (sIn < 1) {
    mprinterr("Error: Bootstrap sample size must be > 0\n");
    return 1;
  }
  if (sIn >= (int)data_->Size()) {
    mprinterr("Error: Bootstrap sample size (%i) must be less than data set size (%zu)\n",
              sIn, data_->Size());
    return 1;
  }
  sample_size_ = sIn;

  if (nresampleIn < 1) {
    mprinterr("Error: Bootstrap number of resamples must be > 0\n");
    return 1;
  }
  if (nresampleIn > (int)data_->Size() - sample_size_)
    mprintf("Warning: Bootstrap # resamples (%i) > data size (%zu) - sample size (%i)\n",
            nresampleIn, data_->Size(), sample_size_);
  n_resample_ = nresampleIn;

  debug_ = debugIn;
  RN_.rn_set(seedIn);
  return 0;
}

// Bootstrap::Resample()
double Bootstrap::Resample(double& Mean) {
  if (data_ == 0) {
    mprinterr("Error: Bootstrap has not been properly initialized.\n");
    return -1.0;
  }

  DataSet_1D const& ds = static_cast<DataSet_1D const&>( *data_ );
  // True if original point already chosen this round
  std::vector<bool> chosen;
  // Hold averages for each resample
  std::vector<double> Avgs(n_resample_, 0.0);
  // Hold average of all resample averages
  Mean = 0.0;
  double d_ndata = (double)ds.Size();
  for (int nsample = 0; nsample != n_resample_; nsample++)
  {
    chosen.assign(ds.Size(), false);
    // Sum up sample_size randomly chosen points
    for (int npoint = 0; npoint != sample_size_; npoint++)
    {
      bool pointOK = false;
      unsigned int pt = 0;
      while (!pointOK)
      {
        pt = (unsigned int)(RN_.rn_gen() * d_ndata);
        pointOK = (chosen[pt] == false);
      }
      chosen[pt] = true;
      Avgs[nsample] += ds.Dval( pt );
    }
    Avgs[nsample] /= (double)sample_size_;
    Mean += Avgs[nsample];
  }

  // Get the standard deviation of all the averages
  Mean /= (double)n_resample_;
  if (debug_ > 0) {
    mprintf("Original mean= %g\n", ds.Avg());
    mprintf("Mean of all resamples= %g\n", Mean);
  }
  double sumdiff2 = 0.0;
  for (int nsample = 0; nsample != n_resample_; nsample++)
  {
    if (debug_ > 1) mprintf("\tMean %i = %g\n", nsample, Avgs[nsample]);
    double diff = Mean - Avgs[nsample];
    sumdiff2 += (diff * diff);
  }
  sumdiff2 /= (double)n_resample_;
  double SD = sqrt(sumdiff2);
  if (debug_ > 0)
    mprintf("SD of all resamples= %g\n", SD);
  return SD;
}
