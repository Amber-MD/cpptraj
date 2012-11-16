// DataSet_float
#include <algorithm> // find
#include "DataSet_float.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet_float::DataSet_float() :
  DataSet(FLOAT, 8, 3, 1)
{}

// DataSet_float::Allocate()
/** Reserve space in the Data and Frames arrays. */
int DataSet_float::Allocate( int sizeIn ) {
  Data_.reserve( sizeIn );
  Frames_.reserve( sizeIn );
  return 0;
}

// DataSet_float::Size()
int DataSet_float::Size() {
  return (int)Data_.size();
}

// DataSet_float::Xmax(()
/** Return the maximum X value added to this set. By convention this is 
  * always the last value added.
  */
int DataSet_float::Xmax() {
  // FIXME: Using this to initialize iterators for buffered write. Should be
  //        a separate routine.
  frame_ = Frames_.begin();
  datum_ = Data_.begin();
  // If no data has been added return 0
  if (Data_.empty()) return 0;
  return Frames_.back();
} 

// DataSet_float::Add()
/** Insert data vIn at frame. */
void DataSet_float::Add(int frame, void *vIn) {
  // Always insert at the end
  // NOTE: No check for duplicate frame values.
  Frames_.push_back( frame );
  Data_.push_back( *((float*)vIn) );
}

// DataSet_float::CurrentDval()
double DataSet_float::CurrentDval() {
  return (double)Data_.back();
}

// DataSet_float::Dval()
double DataSet_float::Dval(int idx) {
  if (idx < 0 || idx >= (int)Data_.size())
    return 0;
  return (double)Data_[idx];
}

// DataSet_float::FrameIsEmpty()
int DataSet_float::FrameIsEmpty(int frame) {
  if ( find( Frames_.begin(), Frames_.end(), frame ) == Frames_.end() )
    return 1;
  return 0;
}

// DataSet_float::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_float::WriteBuffer(CpptrajFile &cbuffer, int frame) {
  while ( frame_ != Frames_.end() && frame > *frame_ ) 
    ++frame_;

  if (frame_ == Frames_.end() || frame != *frame_)
    cbuffer.Printf(data_format_, 0.0);
  else {
    cbuffer.Printf(data_format_, *datum_);
    ++datum_;
    ++frame_;
  }
}

// DataSet_float::Sync()
/** First, non-master threads convert their vectors into C-arrays.
  * These arrays are then sent to the master, where they are put 
  * into the master arrays. It is assumed that master (rank 0) has 
  * first chunk of data, rank 1 has next and so on.
  */
int DataSet_float::Sync() {
  unsigned int masterSize = 0;
  unsigned int dataSize;
  float* values = 0;
  int* frames = 0;

  if (worldsize==1) return 0;

  for ( int rank = 1; rank < worldsize; ++rank) {
    // ----- RANK -------
    if ( worldrank == rank ) {
      // Get size of data on rank.
      dataSize = Data_.size();
      // Send rank size to master
      parallel_sendMaster(&dataSize, 1, rank, 0);
      // If size is 0 on rank, skip this rank.
      if (dataSize == 0) continue;
      // Allocate space on rank
      values = new float[ dataSize ];
      frames = new int[ dataSize ];
      // Send arrays to master
      parallel_sendMaster(frames, dataSize, rank, 0);
      parallel_sendMaster(values, dataSize, rank, 1);
      // Free arrays on rank
      delete[] values;
      delete[] frames;

    // ----- MASTER -----
    } else if (worldrank == 0) {
      // Master receives size from rank
      parallel_sendMaster(&dataSize, 1, rank, 0);
      // If size was 0 on rank, skip rank.
      if (dataSize == 0) continue;
      // Reallocate if necessary
      if (dataSize > masterSize) {
        if ( values != 0 ) delete[] values;
        if ( frames != 0 ) delete[] frames;
        values = new float[ dataSize ];
        frames = new int[ dataSize ];
        masterSize = dataSize;
      }
      // Master receives arrays
      parallel_sendMaster(frames, dataSize, rank, 0);
      parallel_sendMaster(values, dataSize, rank, 1);
      // Insert frames and values to master arrays
      for (unsigned int i = 0; i < dataSize; ++i) {
        Frames_.push_back( frames[i] );
        Data_.push_back( values[i] );
      }
    }
  } // End loop over ranks > 0

  // Free master arrays
  if (worldrank == 0) {
    if ( values != 0 ) delete[] values;
    if ( frames != 0 ) delete[] frames;
  }

  return 0;
}

