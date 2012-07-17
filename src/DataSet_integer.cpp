// DataSet_integer
#include <algorithm> // find
#include "DataSet_integer.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet_integer::DataSet_integer() {
  width_ = 12;
  dType_ = INT;
  SetDataSetFormat(false);
}

// DataSet_integer::Allocate()
/** Reserve space in the Data and Frames arrays. */
int DataSet_integer::Allocate( int sizeIn ) {
  Data_.reserve( sizeIn );
  Frames_.reserve( sizeIn );
  return 0;
}

// DataSet_integer::Size()
int DataSet_integer::Size() {
  return (int)Data_.size();
}

// DataSet_integer::Xmax(()
/** Return the maximum X value added to this set. By convention this is 
  * always the last value added.
  */
int DataSet_integer::Xmax() {
  // FIXME: Using this to initialize iterators for buffered write. Should be
  //        a separate routine.
  frame_ = Frames_.begin();
  datum_ = Data_.begin();
  // If no data has been added return 0
  if (Data_.empty()) return 0;
  return Frames_.back();
} 

// DataSet_integer::Add()
/** Insert data vIn at frame. */
void DataSet_integer::Add(int frame, void *vIn) {
  // Always insert at the end
  // NOTE: No check for duplicate frame values.
  Frames_.push_back( frame );
  Data_.push_back( *((int*)vIn) );
}

// DataSet_integer::CurrentDval()
double DataSet_integer::CurrentDval() {
  return (double)Data_.back();
}

// DataSet_integer::Dval()
double DataSet_integer::Dval(int idx) {
  if (idx < 0 || idx >= (int)Data_.size())
    return 0;
  return (double)Data_[idx];
}

// DataSet_integer::FrameIsEmpty()
int DataSet_integer::FrameIsEmpty(int frame) {
  if ( find( Frames_.begin(), Frames_.end(), frame ) == Frames_.end() )
    return 1;
  return 0;
}

// DataSet_integer::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_integer::WriteBuffer(CharBuffer &cbuffer, int frame) {
  int dval;

  while ( frame_ != Frames_.end() && frame > *frame_ ) 
    ++frame_;

  if (frame_ == Frames_.end() || frame != *frame_)
    dval = 0;
  else {
    dval = *datum_;
    ++datum_;
    ++frame_;
  }

  cbuffer.WriteInteger(data_format_, dval);
}

// DataSet_integer::Width()
int DataSet_integer::Width() {
  return (width_ + leadingSpace_);
}

// DataSet_integer::Sync()
/** First, non-master threads convert their vectors into C-arrays.
  * These arrays are then sent to the master, where they are put 
  * into the master arrays. It is assumed that master (rank 0) has 
  * first chunk of data, rank 1 has next and so on.
  */
int DataSet_integer::Sync() {
  unsigned int masterSize = 0;
  unsigned int dataSize;
  int* values = 0;
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
      values = new int[ dataSize ];
      frames = new int[ dataSize ];
      // Send arrays to master
      parallel_sendMaster(frames, dataSize, rank, 0);
      parallel_sendMaster(values, dataSize, rank, 0);
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
        values = new int[ dataSize ];
        frames = new int[ dataSize ];
        masterSize = dataSize;
      }
      // Master receives arrays
      parallel_sendMaster(frames, dataSize, rank, 0);
      parallel_sendMaster(values, dataSize, rank, 0);
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

