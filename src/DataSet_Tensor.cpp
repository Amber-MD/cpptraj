#include <algorithm> // std::copy
#include "DataSet_Tensor.h"

/** CONSTRUCTOR. */
DataSet_Tensor::DataSet_Tensor() :
  DataSet(TENSOR, GENERIC, TextFormat(TextFormat::DOUBLE, 8, 4, 6), 1)
{}

#ifdef MPI
int DataSet_Tensor::Sync(size_t, std::vector<int> const&, Parallel::Comm const&)
{
  return 1;
}
#endif

/** Reserve memory for given input size. */
int DataSet_Tensor::Allocate( SizeArray const& sizeIn ) {
  if (!sizeIn.empty()) {
    Xvals_.reserve( sizeIn[0] );
    Data_.reserve( sizeIn[0] );
  }
  return 0;
}

/** Add index and tensor. */
void DataSet_Tensor::Add(size_t frame, const void* tIn) {
  Xvals_.push_back( (double)frame );
  Data_.push_back( Ttype( (const double*)tIn ) );
}

/** Write tensor. */
void DataSet_Tensor::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= Data_.size()) {
    cbuffer.Printf(format_.fmt(), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  } else {
    const double* ptr = Data_[ pIn[0] ].Ptr();
    cbuffer.Printf(format_.fmt(), ptr[0], ptr[1], ptr[2],
                                  ptr[3], ptr[4], ptr[5]);
  }
}

/** \return specified index value. */
double DataSet_Tensor::Coord(unsigned int d, size_t p) const {
  if (d > 0) return 0.0;
  return Xvals_[p];
}

/** Append given Tensor set to this one. */
int DataSet_Tensor::Append(DataSet* dsIn) {
  if (dsIn == 0 || dsIn->Empty()) return 0;
  if (dsIn->Type() != TENSOR) return 1;
  Tarray const& tIn = ((DataSet_Tensor*)dsIn)->Data_;
  Farray const& fIn = ((DataSet_Tensor*)dsIn)->Xvals_;
  size_t oldsize = Data_.size();
  Data_.resize( oldsize + tIn.size() );
  std::copy( tIn.begin(), tIn.end(), Data_.begin() + oldsize );
  Xvals_.resize( oldsize + fIn.size() ); // TODO check that tIn and fIn size match?
  std::copy( fIn.begin(), fIn.end(), Xvals_.begin() + oldsize );

  return 0;
}

/** \return Memory usage in bytes. */
size_t DataSet_Tensor::MemUsageInBytes() const {
  return ( sizeof(Tarray) +
           sizeof(Farray) +
           (Data_.size() + Xvals_.size()) * sizeof(double) );
}
