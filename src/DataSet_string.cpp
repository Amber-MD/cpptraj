// DataSet_string
#include "DataSet_string.h"
#ifdef MPI
#include "MpiRoutines.h"
#endif

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
    if (Data_[pIn[0]].size() >= CpptrajFile::BUF_SIZE) {
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

// DataSet_string::Sync()
int DataSet_string::Sync(size_t total, std::vector<int> const& rank_frames) {
# ifdef MPI
  if (worldsize == 1) return 0;
  if (worldrank == 0) {
    // MASTER
    char* block = 0;
    int currentBlockSize = 0;
    size_t pos = Data_.size();
    // Need to increase size of Data on master by number of frames on each other rank.
    int additional_frames = (int)total - rank_frames[0];
    Data_.resize( Data_.size() + additional_frames );
    // Receive data from each rank.
    for (int rank = 1; rank < worldsize; rank++) {
      // Need the size of the char* block
      int blockSize = 0;
      parallel_sendMaster( &blockSize, 1, rank, PARA_INT );
      if (blockSize > currentBlockSize) {
        if (block != 0) delete[] block;
        block = new char[ blockSize ];
        currentBlockSize = blockSize;
      }
      // Receive the block of text
      parallel_sendMaster( block, blockSize, rank, PARA_CHAR );
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
    // Copy each string (including null char) to array
    char* ptr = block;
    for (std::vector<std::string>::const_iterator it = Data_.begin(); it != Data_.end(); ++it)
    {
      size_t len = it->copy( ptr, it->size() );
      ptr[len] = '\0';
      ptr += (len + 1);
    }
    // Send array to master
    parallel_sendMaster( block, blockSize, worldrank, PARA_CHAR );
    delete[] block;
  }
/*
  unsigned int dataSize;
  unsigned int masterStringSize = 0;
  unsigned int stringSize;
  char* values = 0;

  if (worldsize==1) return 0;

  for ( int rank = 1; rank < worldsize; ++rank) {
    if ( worldrank == rank ) {
      // ----- RANK -------
      // Get size of data on rank.
      dataSize = Data_.size();
      // Send rank size to master
      parallel_sendMaster(&dataSize, 1, rank, PARA_INT);
      // If size is 0 on rank, skip this rank.
      if (dataSize == 0) continue;
      // Get sum size of each string on rank (incl. null char).
      stringSize = 0;
      for ( std::vector<std::string>::iterator str_it = Data_.begin(); 
                                          str_it != Data_.end(); ++str_it)
        stringSize += ( (*str_it).size() + 1 ); // +1 for null char.
      // Send sum string size to master
      parallel_sendMaster(&stringSize, 1, rank, PARA_INT);
      // Allocate space on rank
      values = new char[ stringSize ];
      // Copy each string (incl. null char) to the char array
      char* ptr = values;
      for ( std::vector<std::string>::iterator str_it = Data_.begin(); 
                                          str_it != Data_.end(); ++str_it) 
      {
        size_t length = (*str_it).copy( ptr, (*str_it).size() + 1 );
        ptr += length;
      }
      // Send arrays to master
      //parallel_sendMaster(frames, dataSize, rank, PARA_INT);
      parallel_sendMaster(values, stringSize, rank, PARA_CHAR);
      // Free arrays on rank
      delete[] values;
    } else if (worldrank == 0) {
      // ----- MASTER -----
      // Master receives size from rank
      parallel_sendMaster(&dataSize, 1, rank, PARA_INT);
      // If size was 0 on rank, skip rank.
      if (dataSize == 0) continue;
      // Master receives sum string size from rank
      parallel_sendMaster(&stringSize, 1, rank, PARA_INT);
      // Reallocate if necessary
      //if (dataSize > masterSize) {
      //  if ( frames != 0 ) delete[] frames;
      //  frames = new int[ dataSize ];
      //  masterSize = dataSize;
      //}
      if (stringSize > masterStringSize) {
        if ( values != 0 ) delete[] values;
        values = new char[ stringSize ];
        masterStringSize = stringSize;
      }
      // Master receives arrays
      //parallel_sendMaster(frames, dataSize, rank, PARA_INT);
      parallel_sendMaster(values, stringSize, rank, PARA_CHAR);
      // Insert frames and values to master arrays
      char* ptr = values;
      for (unsigned int i = 0; i < dataSize; ++i) {
        Data_.push_back( ptr );
        ptr += ( Data_.back().size() + 1 );
      }
    }
  } // End loop over ranks > 0

  // Free master array
  if (worldrank == 0 && values != 0 ) delete[] values;
*/
# endif
  return 0;
}
