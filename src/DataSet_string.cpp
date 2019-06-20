#include "DataSet_string.h"

size_t DataSet_string::MemUsageInBytes() const {
  size_t mySize = 0;
  for (std::vector<std::string>::const_iterator it = Data_.begin();
                                                it != Data_.end(); ++it)
    mySize += (it->size() * sizeof(char));
  return mySize;
}

// DataSet_string::Allocate()
/** Reserve space in the Data and Frames arrays. */
int DataSet_string::Allocate( SizeArray const& sizeIn ) {
  if (!sizeIn.empty())
    Data_.reserve( sizeIn[0] );
  return 0;
}

// DataSet_string::Add()
/** Insert data vIn at frame. */
void DataSet_string::Add(size_t frame, const void* vIn) {
  if (frame > Data_.size())
    Data_.resize( frame, "NoData" );
  std::string Temp( (const char*)vIn );
  // Check string width. Update format width if necessary.
  if ( (int)Temp.size() > format_.Width() )
    format_.SetWidth(Temp.size());
  // Always insert at the end
  // NOTE: No check for duplicate frame values.
  Data_.push_back( Temp );
}

// DataSet_string::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_string::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= Data_.size())
    cbuffer.Printf(format_.fmt(), "NoData");
  else {
    // Protect against CpptrajFile buffer overflow.
    size_t maxWidth = std::max( (size_t)format_.Width(), Data_[pIn[0]].size() );
    if (maxWidth >= CpptrajFile::BUF_SIZE) {
      // FIXME: Data sets should not have to worry about spaces in format strings.
      if (format_.fmt()[0] == ' ') cbuffer.Printf(" ");
      cbuffer.Write(Data_[pIn[0]].c_str(), Data_[pIn[0]].size());
    } else 
      cbuffer.Printf(format_.fmt(), Data_[pIn[0]].c_str());
  }
}

int DataSet_string::Append(DataSet* dsIn) {
  if (dsIn->Empty()) return 0;
  if (dsIn->Type() != STRING) return 1;
  std::vector<std::string> const& dataIn = 
    static_cast<std::vector<std::string> const&>( ((DataSet_string*)dsIn)->Data() );
  size_t oldsize = Size();
  Data_.resize( oldsize + dataIn.size() );
  std::copy( dataIn.begin(), dataIn.end(), Data_.begin() + oldsize );
  return 0;
}

#ifdef MPI
// DataSet_string::Sync()
int DataSet_string::Sync(size_t total, std::vector<int> const& rank_frames,
                         Parallel::Comm const& commIn)
{
  if (commIn.Size() == 1) return 0;
  if (commIn.Master()) {
    // MASTER
    char* block = 0;
    int currentBlockSize = 0;
    size_t pos = Data_.size();
    // Need to increase size of Data on master by number of frames on each other rank.
    int additional_frames = (int)total - rank_frames[0];
    Data_.resize( Data_.size() + additional_frames );
    // Receive data from each rank.
    for (int rank = 1; rank < commIn.Size(); rank++) {
      // Need the size of the char* block
      int blockSize = 0;
      commIn.SendMaster( &blockSize, 1, rank, MPI_INT );
      if (blockSize > currentBlockSize) {
        if (block != 0) delete[] block;
        block = new char[ blockSize ];
        currentBlockSize = blockSize;
      }
      // Receive the block of text
      commIn.SendMaster( block, blockSize, rank, MPI_CHAR );
      // Convert text block into strings
      const char* ptr = block;
      const char* endptr = block + blockSize;
      while (ptr < endptr) {
        Data_[pos] = std::string(ptr);
        ptr += Data_[pos++].size() + 1; // +1 for null char.
      }
    }
    if (block != 0) delete[] block;
  } else {
    // RANK
    // Get sum size of each string on rank (including null char).
    size_t blockSize = 0;
    for (std::vector<std::string>::const_iterator it = Data_.begin(); it != Data_.end(); ++it)
      blockSize += (it->size() + 1); // +1 for null char.
    char* block = new char[ blockSize ];
    commIn.SendMaster( &blockSize, 1, commIn.Rank(), MPI_INT );
    // Copy each string (including null char) to array
    char* ptr = block;
    for (std::vector<std::string>::const_iterator it = Data_.begin(); it != Data_.end(); ++it)
    {
      size_t len = it->copy( ptr, it->size() );
      ptr[len] = '\0';
      ptr += (len + 1);
    }
    // Send array to master
    commIn.SendMaster( block, blockSize, commIn.Rank(), MPI_CHAR );
    delete[] block;
  }
  return 0;
}
#endif
