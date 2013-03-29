// DataSet_double
#include "DataSet_double.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"

// DataSet_double::Allocate()
/** Reserve space in the Data and Frames arrays. */
int DataSet_double::Allocate1D( size_t sizeIn ) {
  Data_.reserve( sizeIn );
  return 0;
}

// DataSet_double::Add()
/** Insert data vIn at frame. */
void DataSet_double::Add(size_t frame, const void* vIn) {
  if (frame > Data_.size())
    Data_.resize( frame, 0.0 );
  // Always insert at the end
  // NOTE: No check for duplicate frame values.
  Data_.push_back( *((double*)vIn) );
}

// DataSet_double::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_double::WriteBuffer(CpptrajFile &cbuffer, size_t frame) const {
  if (frame >= Data_.size())
    cbuffer.Printf(data_format_, 0.0);
  else
    cbuffer.Printf(data_format_, Data_[frame]);
}

// DataSet_double::Sync()
/** First, non-master threads convert their vectors into C-arrays.
  * These arrays are then sent to the master, where they are put 
  * into the master arrays. It is assumed that master (rank 0) has 
  * first chunk of data, rank 1 has next and so on.
  */
int DataSet_double::Sync() {
  unsigned int dataSize;
  unsigned int masterSize = 0;
  double* values = 0;

  if (worldsize==1) return 0;

  for ( int rank = 1; rank < worldsize; ++rank) {
    if ( worldrank == rank ) {
      // ----- RANK -------
      // Get size of data on rank.
      dataSize = Data_.size();
      // Send rank size to master
      parallel_sendMaster(&dataSize, 1, rank, 0);
      // If size is 0 on rank, skip this rank.
      if (dataSize == 0) continue;
      // Allocate space for temp array on rank, put Data_ into values.
      values = new double[ dataSize ];
      for (size_t i = 0; i < dataSize; ++i)
        values[i] = Data_[i];
      // Send temp array to master
      parallel_sendMaster(values, dataSize, rank, 1);
      // Free array on rank
      delete[] values;
    } else if (worldrank == 0) {
      // ----- MASTER -----
      // Master receives size from rank
      parallel_sendMaster(&dataSize, 1, rank, 0);
      // If size was 0 on rank, skip rank.
      if (dataSize == 0) continue;
      // Reallocate temp array on master if necessary
      if (dataSize > masterSize) {
        if ( values != 0 ) delete[] values;
        values = new double[ dataSize ];
        masterSize = dataSize;
      }
      // Master receives temp array
      parallel_sendMaster(values, dataSize, rank, 1);
      // Insert values to master array
      for (unsigned int i = 0; i < dataSize; ++i)
        Data_.push_back( values[i] );
    }
  } // End loop over ranks > 0

  // Free master array
  if (worldrank == 0 && values != 0 ) delete[] values;

  return 0;
}
