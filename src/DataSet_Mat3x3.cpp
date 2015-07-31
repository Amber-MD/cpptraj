#include "DataSet_Mat3x3.h"

int DataSet_Mat3x3::Allocate(SizeArray const& Nin) {
  if (!Nin.empty())
    data_.reserve( Nin[0] );
  return 0;
}

void DataSet_Mat3x3::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= data_.size()) {
    for (unsigned int i = 0; i < 9; i++)
      cbuffer.Printf(data_format_, 0.0);
  } else {
    for (unsigned int i = 0; i < 9; i++)
      cbuffer.Printf(data_format_, data_[pIn[0]][i]);
  }
}
