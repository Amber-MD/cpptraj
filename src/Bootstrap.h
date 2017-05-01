#ifndef INC_BOOTSTRAP_H
#define INC_BOOTSTRAP_H
#include "Random.h"
#include "DataSet_1D.h"
class Bootstrap {
  public:
    Bootstrap() : sample_size_(0), n_resample_(0), debug_(0) {}
    /// Initialize: DataSet, sample size, number of resamples, random seed, debug
    int Init(DataSet_1D*, int, int, int, int);
    /// \return Error estimation.
    double Resample();
  private:
    Random_Number RN_; ///< Random number generator
    DataSet_1D* data_;
    int sample_size_;
    int n_resample_;
    int debug_;
};
#endif
