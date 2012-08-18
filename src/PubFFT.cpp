#include <cstring> // memset, memcpy
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
  saved_factors_size_(0),
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

void PubFFT::Forward(double* fft_array) {
  cfftf_( fft_size_, fft_array, saved_work_, saved_factors_ );
}

void PubFFT::Back(double* fft_array) {
  cfftb_( fft_size_, fft_array, saved_work_, saved_factors_);
}

