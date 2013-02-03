#include <cstring>
#include "ComplexArray.h"

ComplexArray::ComplexArray(int i) : 
  ndata_(i*2), ncomplex_(i) 
{ 
  if (ndata_ > 0) { 
    data_ = new double[ ndata_ ];
    memset(data_, 0, ndata_ * sizeof(double));
  } else
    data_ = 0;
}

ComplexArray::ComplexArray(ComplexArray const& rhs) :
  ndata_(rhs.ndata_),
  ncomplex_(rhs.ncomplex_)
{
  if (ndata_ > 0) { 
    data_ = new double[ ndata_ ];
    memcpy(data_, rhs.data_, ndata_ * sizeof(double));
  } else
    data_ = 0;
}

ComplexArray& ComplexArray::operator=(ComplexArray const& rhs) {
  if (&rhs == this) return *this;
  if (data_!=0) delete[] data_;
  data_ = 0;
  ndata_ = rhs.ndata_;
  ncomplex_ = rhs.ncomplex_;
  if (ndata_ > 0) {
    data_ = new double[ ndata_ ];
    memcpy(data_, rhs.data_, ndata_ * sizeof(double));
  }
  return *this;
}

ComplexArray::~ComplexArray() { if (data_ != 0) delete[] data_; }

void ComplexArray::Allocate(int i) {
  ndata_ = i * 2;
  ncomplex_ = i;
  if (data_!=0) delete[] data_;
  if (ndata_ > 0)
    data_ = new double[ ndata_ ];
  else
    data_ = 0;
}

void ComplexArray::PadWithZero(int start) {
  int startIdx = 2 * start;
  memset(data_ + startIdx, 0, (ndata_ - startIdx) * sizeof(double));
}

void ComplexArray::Normalize(double norm) {
  for (int i = 0; i < ndata_; ++i)
    data_[i] *= norm;
}

void ComplexArray::SquareModulus() {
  for (int i = 0; i < ndata_; i+=2) {
    data_[i  ] = data_[i  ] * data_[i  ] + data_[i+1] * data_[i+1];
    data_[i+1] = 0.0;
  }
}

// TODO: Check size of rhs
void ComplexArray::ComplexConjTimes(ComplexArray const& rhs) {
  for (int i = 0; i < ndata_; i+=2) {
    double dtmp = data_[i  ] * rhs.data_[i  ] + data_[i+1] * rhs.data_[i+1];
    data_[i+1]  = data_[i  ] * rhs.data_[i+1] - data_[i+1] * rhs.data_[i  ];
    data_[i  ]  = dtmp;
  }
}

// -----------------------------------------------------------------------------
void CorrF_Direct::Allocate(int stepsIn) {
  nsteps_ = stepsIn;
  table_.resize(2*nsteps_, 0.0);
}

void CorrF_Direct::AutoCorr(ComplexArray& data1) {
  int ndata2 = data1.size();
  for (int i = 0; i < ndata2; i++) {
    double dsum = 0.0;
    for (int j = i; j < ndata2; j++) {
      int ind1 = 2 * j;
      int ind2 = 2 * (j-i);
      dsum += data1[ind1] * data1[ind2] + data1[ind1+1] * data1[ind2+1];
    }
    if (i < nsteps_) {
      int ind1 = 2 * i;
      table_[ind1  ] = dsum;
      table_[ind1+1] = 0.0;
    } else
      break;
  }
  std::copy(table_.begin(), table_.end(), data1.CAptr());
}

void CorrF_Direct::CrossCorr(ComplexArray& data1, ComplexArray const& data2) {
  if (data2.size() < data1.size()) return;
  int ndata2 = data1.size();
  for (int i = 0; i < ndata2; i++) {
    double dsum = 0.0;
    double dsumi = 0.0;
    for (int j = i; j < ndata2; j++) {
      int ind1 = 2 * j;
      int ind2 = 2 * (j-i);
      dsum  += data2[ind1] * data1[ind2  ] + data2[ind1+1] * data1[ind2+1];
      dsumi += data2[ind1] * data1[ind2+1] - data2[ind1+1] * data1[ind2  ];
    }
    if(i < nsteps_) {
      int ind1 = 2 * i;
      table_[ind1  ] = dsum;
      table_[ind1+1] = dsumi;
    } else
      break;
  }
  std::copy(table_.begin(), table_.end(), data1.CAptr());
}

// -----------------------------------------------------------------------------
void CorrF_FFT::Allocate(int stepsIn) {
  pubfft_.SetupFFT_NextPowerOf2( stepsIn );
}

void CorrF_FFT::AutoCorr(ComplexArray& data1) {
  pubfft_.Forward( data1.CAptr() );
  // Calculate square modulus of F(data1)
  data1.SquareModulus();
  // Inverse FFT
  pubfft_.Back( data1.CAptr() );
  // Normalize with fft_size (since not done in inverse FFT routine)
  data1.Normalize( 1.0 / (double)pubfft_.size() );
}

void CorrF_FFT::CrossCorr(ComplexArray& data1, ComplexArray& data2) {
  // Cross-correlation
  pubfft_.Forward( data1.CAptr() );
  pubfft_.Forward( data2.CAptr() );
  // Calculate [data1]* x [data2] where * denotes complex conjugate.
  data1.ComplexConjTimes(data2);
  // Inverse FFT
  pubfft_.Back( data1.CAptr() );
  // Normalize with fft_size (since not done in inverse FFT routine)
  data1.Normalize( 1.0 / (double)pubfft_.size() );
}
