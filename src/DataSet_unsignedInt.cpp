#include "DataSet_unsignedInt.h"
//#inc lude "CpptrajStdio.h" // DEBUG

/** Reserve space in the Data and Frames arrays. */
int DataSet_unsignedInt::Allocate( SizeArray const& sizeIn ) {
  if (!sizeIn.empty())
    Data_.reserve( sizeIn[0] );
  return 0;
}

/** Allocate space in Data_ array. */
int DataSet_unsignedInt::MemAlloc( SizeArray const& sizeIn ) {
  if (!sizeIn.empty())
    Data_.resize( sizeIn[0] );
  return 0;
}

/** Copy a block of data of nelts elements from set dptrIn at position pos to startIdx. */
void DataSet_unsignedInt::CopyBlock(size_t startIdx, DataSet const* dptrIn, size_t pos, size_t nelts)
{
  DataSet_unsignedInt const& setIn = static_cast<DataSet_unsignedInt const&>( *dptrIn );
  const unsigned int* ptr = (&(setIn.Data_[0])+pos);
  std::copy( ptr, ptr + nelts, &(Data_[0]) + startIdx );
}

/** Insert data vIn at frame. */
void DataSet_unsignedInt::Add(size_t frame, const void* vIn) {
  if (frame > Data_.size())
    Data_.resize( frame, 0 );
  // Always insert at the end
  // NOTE: No check for duplicate frame values.
  Data_.push_back( *((unsigned int*)vIn) );
}

/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_unsignedInt::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= Data_.size())
    cbuffer.Printf(format_.fmt(), 0);
  else
    cbuffer.Printf(format_.fmt(), Data_[pIn[0]]);
}

/** Append set of same kind. */
int DataSet_unsignedInt::Append(DataSet* dsIn) {
  if (dsIn->Empty()) return 0;
  if (dsIn->Group() != SCALAR_1D) return 1;
  if (dsIn->Type() == UNSIGNED_INTEGER) {
    size_t oldsize = Size();
    std::vector<unsigned int> const& dataIn = ((DataSet_unsignedInt*)dsIn)->Data_;
    Data_.resize( oldsize + dataIn.size() );
    std::copy( dataIn.begin(), dataIn.end(), Data_.begin() + oldsize );
  } else {
    DataSet_1D const& ds = static_cast<DataSet_1D const&>( *dsIn );
    for (unsigned int i = 0; i != ds.Size(); i++)
      Data_.push_back( (unsigned int)ds.Dval(i) );
  }
  return 0;
}

#ifdef MPI
// DataSet_unsignedInt::Sync()
int DataSet_unsignedInt::Sync(size_t total, std::vector<int> const& rank_frames,
                          Parallel::Comm const& commIn)
{
  if (commIn.Size()==1) return 0;
  if (commIn.Master()) {
    // Resize for total number of frames.
    Data_.resize( total );
    unsigned int* endptr = &(Data_[0]) + rank_frames[0];
    // Receive data from each rank
    for (int rank = 1; rank < commIn.Size(); rank++) {
      commIn.SendMaster( endptr, rank_frames[rank], rank, MPI_UNSIGNED );
      endptr += rank_frames[rank];
    }
  } else // Send data to master //TODO adjust for repeated additions?
    commIn.SendMaster( &(Data_[0]), Data_.size(), commIn.Rank(), MPI_UNSIGNED );
  return 0;
}

/** Receive unsigned integer data from rank, append to after offset. */
int DataSet_unsignedInt::Recv(size_t total, unsigned int offset, int nframes,
                              int fromRank, int tag, Parallel::Comm const& commIn)
{
  Data_.resize( total );
  unsigned int* d_beg = &(Data_[0]) + offset;
  //rprintf("DEBUG: Receive %i frames fromRank %i tag %i\n", nframes, fromRank, tag);
  if (commIn.Recv( d_beg, nframes, MPI_UNSIGNED, fromRank, tag )) return 1;
  // TODO Do SetNeedsSync false here?
  return 0;
}

/** Send unsigned integer data to rank. */
int DataSet_unsignedInt::Send(int toRank, int tag, Parallel::Comm const& commIn)
const
{
  //rprintf("DEBUG: Send %zu frames toRank %i tag %i\n", Data_.size(), toRank, tag);
  return commIn.Send( (void*)&(Data_[0]), Data_.size(), MPI_UNSIGNED, toRank, tag );
}
#endif
