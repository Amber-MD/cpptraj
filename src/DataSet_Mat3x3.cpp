#include "DataSet_Mat3x3.h"

int DataSet_Mat3x3::Allocate1D(size_t Nin) {
  data_.reserve( Nin );
  return 0;
}

void DataSet_Mat3x3::WriteBuffer(CpptrajFile &cbuffer, size_t frameIn) const {
  if (frameIn >= data_.size()) {
    for (unsigned int i = 0; i < 9; i++)
      cbuffer.Printf(data_format_, 0.0);
  } else {
    for (unsigned int i = 0; i < 9; i++)
      cbuffer.Printf(data_format_, data_[frameIn][i]);
  }
}
