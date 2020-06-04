#include "DataSet_Vector_Scalar.h"
#include "CpptrajStdio.h"
#include <algorithm> // std::copy

/// CONSTRUCTOR TODO should this be in VECTOR_1D group?
DataSet_Vector_Scalar::DataSet_Vector_Scalar() :
  DataSet(VECTOR_SCALAR, GENERIC, TextFormat(TextFormat::DOUBLE, 8, 4, 4), 1)
{ }

/** Reserve space in vector array. */
int DataSet_Vector_Scalar::Allocate(SizeArray const& Nin) {
  if (!Nin.empty()) {
    vecs_.reserve( Nin[0] );
    vals_.reserve( Nin[0] );
  }
  return 0;
}

/** Write vector + scalar to output file. */
void DataSet_Vector_Scalar::WriteBuffer(CpptrajFile& cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= Size()) {
    cbuffer.Printf(format_.fmt(), 0.0, 0.0, 0.0, 0.0); /// VXYZ SCALAR
  } else {
    Vec3 const& vxyz = vecs_[pIn[0]];
    cbuffer.Printf(format_.fmt(), vxyz[0], vxyz[1], vxyz[2], vals_[pIn[0]]);
  }
}

/** Append vector scalar set to this set. */
int DataSet_Vector_Scalar::Append(DataSet* dsIn) {
  if (dsIn == 0) return 1;
  if (dsIn->Empty()) return 0;
  if (dsIn->Type() != VECTOR_SCALAR) return 1;
  Varray const& vIn = ((DataSet_Vector_Scalar*)dsIn)->vecs_;
  Darray const& dIn = ((DataSet_Vector_Scalar*)dsIn)->vals_;
  size_t oldsize = vecs_.size();
  vecs_.resize( oldsize + vIn.size() );
  std::copy( vIn.begin(), vIn.end(), vecs_.begin() + oldsize );
  vals_.resize( oldsize + dIn.size() );
  std::copy( dIn.begin(), dIn.end(), vals_.begin() + oldsize );
  return 0;
}

/** \return memory usage of the set in bytes. */
size_t DataSet_Vector_Scalar::MemUsageInBytes() const {
  return ((vecs_.size()*Vec3::DataSize()) + vals_.size());
}
