#include "DataSet_Vector_Scalar.h"
#include "CpptrajStdio.h"
#include <algorithm> // std::copy

/// CONSTRUCTOR TODO should this be in VECTOR_1D group?
DataSet_Vector_Scalar::DataSet_Vector_Scalar() :
  // TODO marking this as 0-dimension so that it is forced to use
  //      the DataIO_Peaks format, but may want to change this down the
  //      road.
  DataSet(VECTOR_SCALAR, GENERIC, TextFormat(TextFormat::DOUBLE, 8, 4, 4), 0)
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

#ifdef MPI
int DataSet_Vector_Scalar::Sync(size_t total, std::vector<int> const& rank_frames,
                         Parallel::Comm const& commIn)
{
  if (commIn.Size()==1) return 0;
  double buf[4];
  if (commIn.Master()) {
    // Resize to accept data from other ranks.
    vecs_.resize( total );
    vals_.resize( total );
    int vidx = rank_frames[0];
    for (int rank = 1; rank < commIn.Size(); rank++) {
      for (int ridx = 0; ridx != rank_frames[rank]; ridx++, vidx++) {
        commIn.SendMaster( buf, 4, rank, MPI_DOUBLE );
        std::copy( buf, buf+3, vecs_[vidx].Dptr() );
        vals_[vidx] = buf[3];
      }
    }
  } else {
    // Send data to master
    for (unsigned int ridx = 0; ridx != vecs_.size(); ++ridx) {
      std::copy( vecs_[ridx].Dptr(), vecs_[ridx].Dptr()+3, buf );
      buf[3] = vals_[ridx];
      commIn.SendMaster( buf, 4, commIn.Rank(), MPI_DOUBLE );
    }
  }
  return 0;
}
#endif
