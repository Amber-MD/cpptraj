#include "DataSet_integer.h"
#ifdef MPI
# include "Parallel.h"
#endif

// DataSet_integer::Allocate()
/** Reserve space in the Data and Frames arrays. */
int DataSet_integer::Allocate( SizeArray const& sizeIn ) {
  if (!sizeIn.empty())
    Data_.reserve( sizeIn[0] );
  return 0;
}

// DataSet_integer::Add()
/** Insert data vIn at frame. */
void DataSet_integer::Add(size_t frame, const void* vIn) {
  if (frame > Data_.size())
    Data_.resize( frame, 0 );
  // Always insert at the end
  // NOTE: No check for duplicate frame values.
  Data_.push_back( *((int*)vIn) );
}

// DataSet_integer::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_integer::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= Data_.size())
    cbuffer.Printf(format_.fmt(), 0);
  else
    cbuffer.Printf(format_.fmt(), Data_[pIn[0]]);
}

int DataSet_integer::Append(DataSet* dsIn) {
  if (dsIn->Empty()) return 0;
  if (dsIn->Group() != SCALAR_1D) return 1;
  if (dsIn->Type() == INTEGER) {
    size_t oldsize = Size();
    std::vector<int> const& dataIn = ((DataSet_integer*)dsIn)->Data_;
    Data_.resize( oldsize + dataIn.size() );
    std::copy( dataIn.begin(), dataIn.end(), Data_.begin() + oldsize );
  } else {
    DataSet_1D const& ds = static_cast<DataSet_1D const&>( *dsIn );
    for (unsigned int i = 0; i != ds.Size(); i++)
      Data_.push_back( (int)ds.Dval(i) );
  }
  return 0;
}

// DataSet_integer::Sync()
int DataSet_integer::Sync(size_t total, std::vector<int> const& rank_frames) {
# ifdef MPI
  if (Parallel::World().Size()==1) return 0;
  if (Parallel::World().Master()) {
    // Resize for total number of frames.
    Data_.resize( total );
    int* endptr = &(Data_[0]) + rank_frames[0];
    // Receive data from each rank
    for (int rank = 1; rank < Parallel::World().Size(); rank++) {
      Parallel::World().SendMaster( endptr, rank_frames[rank], rank, MPI_INT );
      endptr += rank_frames[rank];
    }
  } else // Send data to master //TODO adjust for repeated additions?
    Parallel::World().SendMaster( &(Data_[0]), Data_.size(), Parallel::World().Rank(), MPI_INT );
# endif
  return 0;
}
