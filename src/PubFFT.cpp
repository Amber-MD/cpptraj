#include <cstring> // memset, memcpy
#include <cmath> // log, ldexp
#include "PubFFT.h"

extern "C" {
  // pubfft.F90 functions
  void cffti_(int&, double*, int*);
  void cfftf_(int&, double*, double*, int*);
  void cfftb_(int&, double*, double*, int*);
}

// CONSTRUCTOR
PubFFT::PubFFT() :
  fft_size_(0),
  saved_factors_size_(30),
  saved_work_size_(0),
  saved_factors_(0),
  saved_work_(0)
{}

// DESTRUCTOR
PubFFT::~PubFFT() {
  if (saved_factors_ != 0) delete[] saved_factors_;
  if (saved_work_ != 0) delete[] saved_work_;
}

// CONSTRUCTOR
/// Takes FFT size as input
PubFFT::PubFFT(int fft_sizeIn) :
  fft_size_(fft_sizeIn),
  saved_factors_size_(30),
  saved_work_size_(4 * fft_size_)
{
  saved_factors_ = new int[ saved_factors_size_ ];
  memset(saved_factors_, 0, saved_factors_size_ * sizeof(int));
  saved_work_ = new double[ saved_work_size_ ];
  memset(saved_work_, 0, saved_work_size_ * sizeof(double));
  cffti_( fft_size_, saved_work_, saved_factors_ );
}

// COPY CONSTRUCTOR
PubFFT::PubFFT(const PubFFT& rhs) :
  fft_size_(rhs.fft_size_),
  saved_factors_size_(rhs.saved_factors_size_),
  saved_work_size_(rhs.saved_work_size_),
  saved_factors_(0),
  saved_work_(0)
{
  if (saved_factors_size_ > 0) {
    saved_factors_ = new int[ saved_factors_size_ ];
    memcpy(saved_factors_, rhs.saved_factors_, saved_factors_size_ * sizeof(int));
  }
  if (saved_work_size_ > 0) {
    saved_work_ = new double[ saved_work_size_ ];
    memcpy(saved_work_, rhs.saved_work_, saved_work_size_ * sizeof(double));
  }
}

// Assignment operator
PubFFT& PubFFT::operator=(const PubFFT& rhs) {
  if (this == &rhs) return *this;
  delete[] saved_work_;
  delete[] saved_factors_;
  fft_size_ = rhs.fft_size_;
  saved_factors_size_ = rhs.saved_factors_size_;
  saved_work_size_ = rhs.saved_work_size_;
  if (saved_factors_size_ > 0) {
    saved_factors_ = new int[ saved_factors_size_ ];
    memcpy(saved_factors_, rhs.saved_factors_, saved_factors_size_ * sizeof(int));
  } else
    saved_factors_ = 0;
  if (saved_work_size_ > 0) {
    saved_work_ = new double[ saved_work_size_ ];
    memcpy(saved_work_, rhs.saved_work_, saved_work_size_ * sizeof(double));
  } else
    saved_work_ = 0;
  return *this;
}

// PubFFT::Forward()
void PubFFT::Forward(double* fft_array) {
  cfftf_( fft_size_, fft_array, saved_work_, saved_factors_ );
}

// PubFFT::Back()
void PubFFT::Back(double* fft_array) {
  cfftb_( fft_size_, fft_array, saved_work_, saved_factors_);
}

// PubFFT::SetupFFT()
int PubFFT::SetupFFT(int sizeIn) {
  int ndata = ldexp( 1.0, (int)(log((double)(4 * sizeIn - 1)) / log(2.0)) + 1);
  fft_size_ = ndata / 2;
  saved_work_size_ = 2 * ndata;
  // NOTE: Do not change saved_factors_size
  if (saved_work_ != 0) delete[] saved_work_;
  if (saved_factors_ != 0) delete[] saved_factors_;
  // Reallocate FFT workspace
  if (saved_factors_size_ > 0) {
    saved_factors_ = new int[ saved_factors_size_ ];
    memset(saved_factors_, 0, saved_factors_size_ * sizeof(int));
  } else
    saved_factors_ = 0;
  if (saved_work_size_ > 0) {
    saved_work_ = new double[ saved_work_size_ ];
    memset(saved_work_, 0, saved_work_size_ * sizeof(double));
  } else
    saved_work_ = 0;
  // NOTE: Should this be called if fft_size is 0?
  cffti_( fft_size_, saved_work_, saved_factors_ );
  return 0;
}

// PubFFT::CorF_FFT()
/** /author Original Authors: Alrun N. Koller & H. Gohlke
  *
  * Calculates correlation functions using the Wiener-Khinchin-Theorem
  * (s. Comp. Sim. of Liquids, p. 188)
  *
  * - ndata is the length of the arrays data1/2
  * - data1/2 is a real array of complex numbers: 
  *   data1/2(1)=real1/2(1), data1/2(2)=img1/2(1), ...
  * - it is recommended that the discrete data is appended
  *     with as much zeros to avoid spurious correlations
  * - in addition, data MUST have the dimension of power of 2
  *     (pad the "real" data (plus zeros) with additional 
  *     zeros up to the next power of 2)
  * - the result is not yet normalized by (no_of_discrete_data - t)**-1 (!)
  */
void PubFFT::CorF_FFT(int ndata, double* data1, double* data2) {
  //int ndata = fft_size_ * 2;
  // FFT data
  cfftf_(fft_size_, data1, saved_work_, saved_factors_);
  if(data2 != NULL)
    cfftf_(fft_size_, data2, saved_work_, saved_factors_);

  // Calc square modulus (in case of cross-correlation: calc [F(data1)]' * F(data2),
  //   where [F(data1)]' is the complex conjugate of F(data1)).
  if (data2 == NULL) {
    for (int i = 0; i < ndata; i+=2) {
      data1[i  ] = data1[i  ] * data1[i  ] + data1[i+1] * data1[i+1];
      data1[i+1] = 0.0;
    }
  } else {
    for(int i = 0; i < ndata; i+=2) {
      double dtmp       = data1[i  ] * data2[i  ] + data1[i+1] * data2[i+1];
      data1[i+1] = data1[i  ] * data2[i+1] - data2[i  ] * data1[i+1];
      data1[i  ] = dtmp;
    }
  }

  // Inverse FFT
  cfftb_(fft_size_, data1, saved_work_, saved_factors_);

  // Normalize with fft_size (since not done in inverse FFT routine)
  double norm = 1.0 / ((double) (fft_size_));
  for(int i = 0; i < ndata; ++i)
    data1[i] *= norm;
}
