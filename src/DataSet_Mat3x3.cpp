#include "DataSet_Mat3x3.h"

int DataSet_Mat3x3::Allocate(SizeArray const& Nin) {
  if (!Nin.empty())
    data_.reserve( Nin[0] );
  return 0;
}

void DataSet_Mat3x3::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= data_.size())
    cbuffer.Printf(format_.fmt(), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  else {
    Matrix_3x3 const& m = data_[pIn[0]];
    cbuffer.Printf(format_.fmt(), m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]);
  }
}

int DataSet_Mat3x3::Append(DataSet* dsIn) {
  if (dsIn->Empty()) return 0;
  if (dsIn->Type() != MAT3X3) return 1;
  Marray const& mIn = ((DataSet_Mat3x3*)dsIn)->data_;
  size_t oldsize = Size();
  data_.resize( oldsize + mIn.size() );
  std::copy( mIn.begin(), mIn.end(), data_.begin() + oldsize );
  return 0;
}

#ifdef MPI
int DataSet_Mat3x3::Sync(size_t total, std::vector<int> const& rank_frames,
                         Parallel::Comm const& commIn)
{
  if (commIn.Size()==1) return 0;
  if (commIn.Master()) {
    // Resize to accept data from other ranks.
    data_.resize( total );
    int midx = rank_frames[0]; // Index on master
    for (int rank = 1; rank < commIn.Size(); rank++) {
      for (int ridx = 0; ridx != rank_frames[rank]; ridx++, midx++)
        // TODO: Consolidate to 1 send/recv via arrays?
        commIn.SendMaster( data_[midx].Dptr(), 9, rank, MPI_DOUBLE );
    }
  } else { // Send data to master
    for (unsigned int ridx = 0; ridx != data_.size(); ++ridx)
      commIn.SendMaster( data_[ridx].Dptr(), 9, commIn.Rank(), MPI_DOUBLE );
  }
  return 0;
}
#endif
