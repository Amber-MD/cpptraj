// DataSet_integer
#include "DataSet_integer.h"
#ifdef MPI
#include "MpiRoutines.h"
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
#ifdef MPI
  if (worldsize==1) return 0;
  if (worldrank == 0) {
    int* endptr = &(Data_[0]) + Data_.size();
    // Need to increase size of Data on master by number of frames on each other rank.
    int additional_frames = (int)total - rank_frames[0];
    Data_.resize( Data_.size() + additional_frames );
    // Receive data from each rank.
    for (int rank = 1; rank < worldsize; rank++) {
      parallel_sendMaster( endptr, rank_frames[rank], rank, PARA_INT );
      endptr += rank_frames[rank];
    }
  } else // Send data to master //TODO adjust for repeated additions?
    parallel_sendMaster( &(Data_[0]), Data_.size(), worldrank, PARA_INT );
/*
  unsigned int dataSize;
  unsigned int masterSize = 0;
  int* values = 0;

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
      // Allocate space for temp array on rank, put Data_ into values.
      values = new int[ dataSize ];
      std::copy(Data_.begin(), Data_.end(), values);
      //frames = new int[ dataSize ];
      // Send arrays to master
      //parallel_sendMaster(frames, dataSize, rank, PARA_INT);
      parallel_sendMaster(values, dataSize, rank, PARA_INT);
      // Free arrays on rank
      delete[] values;
    } else if (worldrank == 0) {
      // ----- MASTER -----
      // Master receives size from rank
      parallel_sendMaster(&dataSize, 1, rank, PARA_INT);
      // If size was 0 on rank, skip rank.
      if (dataSize == 0) continue;
      // Reallocate temp array on master if necessary
      if (dataSize > masterSize) {
        if ( values != 0 ) delete[] values;
        values = new int[ dataSize ];
        masterSize = dataSize;
      }
      // Master receives arrays
      //parallel_sendMaster(frames, dataSize, rank, PARA_INT);
      parallel_sendMaster(values, dataSize, rank, PARA_INT);
      // Insert frames and values to master arrays
      for (unsigned int i = 0; i < dataSize; ++i) {
        //Frames_.push_back( frames[i] );
        Data_.push_back( values[i] );
      }
    }
  } // End loop over ranks > 0

  // Free master array
  if (worldrank == 0 && values != 0 ) delete[] values;
*/
#endif
  return 0;
}
