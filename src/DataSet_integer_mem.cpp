#include "DataSet_integer_mem.h"
//#inc lude "CpptrajStdio.h" // DEBUG

// DataSet_integer_mem::Allocate()
/** Reserve space in the Data and Frames arrays. */
int DataSet_integer_mem::Allocate( SizeArray const& sizeIn ) {
  if (!sizeIn.empty())
    Data_.reserve( sizeIn[0] );
  return 0;
}

// DataSet_integer_mem::Add()
/** Insert data vIn at frame. */
void DataSet_integer_mem::Add(size_t frame, const void* vIn) {
  if (frame > Data_.size())
    Data_.resize( frame, 0 );
  // Always insert at the end
  // NOTE: No check for duplicate frame values.
  Data_.push_back( *((int*)vIn) );
}

// DataSet_integer_mem::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_integer_mem::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= Data_.size())
    cbuffer.Printf(format_.fmt(), 0);
  else
    cbuffer.Printf(format_.fmt(), Data_[pIn[0]]);
}

int DataSet_integer_mem::Append(DataSet* dsIn) {
  if (dsIn->Empty()) return 0;
  if (dsIn->Group() != SCALAR_1D) return 1;
  if (dsIn->Type() == INTEGER) { // TODO check disk cache
    size_t oldsize = Size();
    std::vector<int> const& dataIn = ((DataSet_integer_mem*)dsIn)->Data_;
    Data_.resize( oldsize + dataIn.size() );
    std::copy( dataIn.begin(), dataIn.end(), Data_.begin() + oldsize );
  } else {
    DataSet_1D const& ds = static_cast<DataSet_1D const&>( *dsIn );
    for (unsigned int i = 0; i != ds.Size(); i++)
      Data_.push_back( (int)ds.Dval(i) );
  }
  return 0;
}

#ifdef MPI
// DataSet_integer_mem::Sync()
int DataSet_integer_mem::Sync(size_t total, std::vector<int> const& rank_frames,
                          Parallel::Comm const& commIn)
{
  if (commIn.Size()==1) return 0;
  if (commIn.Master()) {
    // Resize for total number of frames.
    Data_.resize( total );
    int* endptr = &(Data_[0]) + rank_frames[0];
    // Receive data from each rank
    for (int rank = 1; rank < commIn.Size(); rank++) {
      commIn.SendMaster( endptr, rank_frames[rank], rank, MPI_INT );
      endptr += rank_frames[rank];
    }
  } else // Send data to master //TODO adjust for repeated additions?
    commIn.SendMaster( &(Data_[0]), Data_.size(), commIn.Rank(), MPI_INT );
  return 0;
}

/** Receive integer data from rank, append to after offset. */
int DataSet_integer_mem::Recv(size_t total, unsigned int offset, int nframes,
                              int fromRank, int tag, Parallel::Comm const& commIn)
{
  Data_.resize( total );
  int* d_beg = &(Data_[0]) + offset;
  //rprintf("DEBUG: Receive %i frames fromRank %i tag %i\n", nframes, fromRank, tag);
  if (commIn.Recv( d_beg, nframes, MPI_INT, fromRank, tag )) return 1;
  // TODO Do SetNeedsSync false here?
  return 0;
}

/** Send integer data to rank. */
int DataSet_integer_mem::Send(int toRank, int tag, Parallel::Comm const& commIn)
const
{
  //rprintf("DEBUG: Send %zu frames toRank %i tag %i\n", Data_.size(), toRank, tag);
  return commIn.Send( (void*)&(Data_[0]), Data_.size(), MPI_INT, toRank, tag );
}
#endif
