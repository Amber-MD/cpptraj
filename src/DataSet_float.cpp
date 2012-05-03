// DataSet_float
#include "DataSet_float.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"
using namespace std;

// CONSTRUCTOR
DataSet_float::DataSet_float() {
  width_ = 8;
  precision_ = 3;
  dType_ = FLOAT;
  SetDataSetFormat(false);
}

// DataSet_float::Size()
int DataSet_float::Size() {
  return (int)Data_.size();
}

// DataSet_float::Begin()
void DataSet_float::Begin() {
  datum_ = Data_.begin();
}

// DataSet_float::NextValue();
bool DataSet_float::NextValue() {
  ++datum_;
  if (datum_ == Data_.end()) return false;
  return true;
}

// DataSet_float::CurrentValue()
double DataSet_float::CurrentValue() {
  return (double)(*datum_).second;
}


// DataSet_float::Xmax(()
/** Return the maximum X value added to this set. By convention this is 
  * always the last value added.
  */
int DataSet_float::Xmax() {
  // If no data has been added return 0
  //if (current_==0) return 0;
  if (Data_.empty()) return 0;
  datum_ = Data_.end();
  --datum_;
  return ( (*datum_).first );
} 

// DataSet_float::Add()
/** Insert data vIn at frame. */
void DataSet_float::Add(int frame, void *vIn) {
  float *value = (float*) vIn;

  // Always insert at the end
  //datum=Data.end();
  //Data.insert( datum, pair<int,float>(frame, *value) );
  Data_[frame] = (*value);
  //++current_;
}

// DataSet_float::Get()
/** Get data at frame, put into vOut. Return 1 if no data at frame.
  */
int DataSet_float::Get(void *vOut, int frame) {
  float *value;
  
  if (vOut==NULL) return 1;
  //mprintf("DEBUG: Attempting to get float frame %i\n",frame);
  value = (float*) vOut;
  datum_ = Data_.find( frame );
  if (datum_ == Data_.end()) return 1;
  //mprintf("DEBUG: Float frame %i is %lf\n",frame,(*datum).second);
  *value = (*datum_).second;
  return 0;
}

// DataSet_float::Dval()
double DataSet_float::Dval(int idx) {
  float val;
  if (Get(&val,idx)) return 0;
  return (double)val;
}

// DataSet_float::FrameIsEmpty()
int DataSet_float::FrameIsEmpty(int frame) {
  datum_ = Data_.find( frame );
  if (datum_ == Data_.end()) return 1;
  return 0;
}

// DataSet_float::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_float::WriteBuffer(CharBuffer &cbuffer, int frame) {
  double dval;
  datum_ = Data_.find( frame );
  if (datum_ == Data_.end())
    dval = 0.0;
  else
    dval = (double)(*datum_).second;
  cbuffer.WriteDouble(data_format_, dval);
}

// DataSet_float::Width()
int DataSet_float::Width() {
  return (width_ + leadingSpace_);
}

// DataSet_float::Sync()
/** Since it seems to be very difficult (or impossible) to define Classes
  * as MPI datatypes, first non-master threads need to convert their maps
  * into 2 arrays, an int array containing frame #s and a double array
  * containing mapped values. These arrays are then sent to the master,
  * where they are converted pairs and inserted into the master map.
  */
int DataSet_float::Sync() {
  int rank, i, dataSize;
  int *Frames;
  float *Values;

  if (worldsize==1) return 0;

  for (rank = 1; rank < worldsize; rank ++) {
    // Get size of map on rank 
    if (worldrank>0) {
      // NOTE: current should be equal to size(). Check for now
      //rprintf("DataSet_float syncing. current=%i, size=%u\n",
      //        current_, Data_.size());
      /*if (current != (int) Data.size()) {
        rprintf("ERROR: current and map size are not equal.\n");
        return 1;
      }*/
      dataSize = (int)Data_.size();
    }

    // Send size of map on rank to master, allocate arrays on rank and master
    parallel_sendMaster(&dataSize, 1, rank, 0);
    rprintf("DataSet_float allocating %i for send/recv\n",dataSize);
    Frames = new int[ dataSize ];
    Values = new float[ dataSize ];
      
    // On non-master convert map to int and double arrays.
    if (worldrank > 0) {
      i=0;
      for ( datum_ = Data_.begin(); datum_ != Data_.end(); datum_++ ) {
        Frames[i]=(*datum_).first;
        Values[i]=(*datum_).second;
        ++i;
      }
    }

    // Send arrays to master
    parallel_sendMaster(Frames, dataSize, rank, 0);
    parallel_sendMaster(Values, dataSize, rank, 1);

    // On master convert arrays to pairs and insert to master map
    if (worldrank==0) {
      for (i=0; i < dataSize; i++) 
        Data_.insert( pair<int,float>( Frames[i], Values[i] ) );
    }

    // Free arrays
    delete[] Frames;
    delete[] Values;
  } // End loop over ranks>0

  return 0;
}

