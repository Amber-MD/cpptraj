// DataSet_string
#include <cstring> // strcpy
#include <algorithm>
#include "DataSet_string.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet_string::DataSet_string() :
  DataSet(STRING, 1, 0, 1)
{}

// DataSet_string::Allocate()
/** Reserve space in the Data and Frames arrays. */
int DataSet_string::Allocate( int sizeIn ) {
  Data_.reserve( sizeIn );
  Frames_.reserve( sizeIn );
  return 0;
}

// DataSet_string::Size()
int DataSet_string::Size() {
  return (int)Data_.size();
}

// DataSet_string::Xmax(()
/** Return the maximum X value added to this set. By convention this is 
  * always the last value added.
  */
int DataSet_string::Xmax() {
  // FIXME: Using this to initialize iterators for buffered write. Should be
  //        a separate routine.
  frame_ = Frames_.begin();
  datum_ = Data_.begin();
  // If no data has been added return 0
  if (Data_.empty()) return 0;
  return Frames_.back();
}

// DataSet_string::Add()
/** Insert data vIn at frame. If the size of the input string is greater
  * than the current width, increase the width.
  * String expects char*
  */
void DataSet_string::Add(int frame, void *vIn) {
  char* value = (char*)vIn;
  std::string Temp( value );
  int strsize = (int)Temp.size();
  if (strsize > width_) width_ = strsize;
  // Always insert at the end
  // NOTE: No check for duplicate frame values.
  Frames_.push_back( frame );
  Data_.push_back( Temp );
}

// DataSet_string::FrameIsEmpty()
int DataSet_string::FrameIsEmpty(int frame) {
  if ( find( Frames_.begin(), Frames_.end(), frame ) == Frames_.end() )
    return 1;
  return 0;
}

// DataSet_string::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write NoData..
  */
void DataSet_string::WriteBuffer(CpptrajFile &cbuffer, int frame) {
  while ( frame_ != Frames_.end() && frame > *frame_ )
    ++frame_;

  if (frame_ == Frames_.end() || frame != *frame_)
    cbuffer.Printf(data_format_,"NoData");
  else {
    cbuffer.Printf(data_format_, (*datum_).c_str());
    ++datum_;
    ++frame_;
  }
}

// DataSet_string::Sync()
/** First, non-master threads convert their vectors into C-arrays.
  * These arrays are then sent to the master, where they are put 
  * into the master arrays. It is assumed that master (rank 0) has 
  * first chunk of data, rank 1 has next and so on.
  */

int DataSet_string::Sync() {
  unsigned int masterSize = 0;
  unsigned int dataSize;
  unsigned int masterStringSize = 0;
  unsigned int stringSize;
  char* values = 0;
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
      // Get sum size of each string on rank (incl. null char).
      stringSize = 0;
      for ( DType::iterator str_it = Data_.begin(); str_it != Data_.end(); ++str_it)
        stringSize += ( (*str_it).size() + 1 ); // +1 for null char.
      // Send sum string size to master
      parallel_sendMaster(&stringSize, 1, rank, 0);
      // Allocate space on rank
      values = new char[ stringSize ];
      frames = new int[ dataSize ];
      // Copy each string (incl. null char) to the char array
      char* ptr = values;
      for ( DType::iterator str_it = Data_.begin(); str_it != Data_.end(); ++str_it) {
        strcpy( ptr, (*str_it).c_str() );
        ptr += ( (*str_it).size() + 1 );
      }
      // Send arrays to master
      parallel_sendMaster(frames, dataSize, rank, 0);
      parallel_sendMaster(values, stringSize, rank, 2);
      // Free arrays on rank
      delete[] values;
      delete[] frames;

    // ----- MASTER -----
    } else if (worldrank == 0) {
      // Master receives size from rank
      parallel_sendMaster(&dataSize, 1, rank, 0);
      // If size was 0 on rank, skip rank.
      if (dataSize == 0) continue;
      // Master receives sum string size from rank
      parallel_sendMaster(&stringSize, 1, rank, 0);
      // Reallocate if necessary
      if (dataSize > masterSize) {
        if ( frames != 0 ) delete[] frames;
        frames = new int[ dataSize ];
        masterSize = dataSize;
      }
      if (stringSize > masterStringSize) {
        if ( values != 0 ) delete[] values;
        values = new char[ stringSize ];
        masterStringSize = stringSize;
      }
      // Master receives arrays
      parallel_sendMaster(frames, dataSize, rank, 0);
      parallel_sendMaster(values, stringSize, rank, 2);
      // Insert frames and values to master arrays
      char* ptr = values;
      for (unsigned int i = 0; i < dataSize; ++i) {
        Frames_.push_back( frames[i] );
        Data_.push_back( ptr );
        ptr += ( Data_.back().size() + 1 );
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

