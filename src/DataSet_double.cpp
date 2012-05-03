// DataSet_double
#include "DataSet_double.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"
using namespace std;

// CONSTRUCTOR
DataSet_double::DataSet_double() {
  width_ = 12;
  precision_ = 4;
  dType_ = DOUBLE;
  SetDataSetFormat(false);
}

// DataSet_double::Size()
int DataSet_double::Size() {
  return (int)Data_.size();
}

// DataSet_double::Begin()
void DataSet_double::Begin() {
  datum_ = Data_.begin();
}

// DataSet_double::NextValue();
bool DataSet_double::NextValue() {
  ++datum_;
  if (datum_ == Data_.end()) return false;
  return true;
}

// DataSet_double::CurrentValue()
double DataSet_double::CurrentValue() {
  return (*datum_).second;
}

// DataSet_double::Xmax(()
/** Return the maximum X value added to this set. By convention this is 
  * always the last value added.
  */
int DataSet_double::Xmax() {
  // If no data has been added return 0
  //if (current_==0) return 0;
  if (Data_.empty()) return 0;
  datum_ = Data_.end();
  --datum_;
  return ( (*datum_).first );
} 

// DataSet_double::Add()
/** Insert data vIn at frame.  */
void DataSet_double::Add(int frame, void *vIn) {
  double *value;

  value = (double*) vIn;

  // Always insert at the end
  Data_[frame] = (*value);
  //it=Data.end();
  //Data.insert( it, pair<int,double>(frame, *value) );
  //current_++;
}

// DataSet_double::Get()
/** Get data at frame, put into vOut. 
  * \return 1 if no data at frame.
  */
int DataSet_double::Get(void *vOut, int frame) {
  double *value;
  
  if (vOut==NULL) return 1;
  //mprintf("DEBUG: Attempting to get double frame %i\n",frame);
  value = (double*) vOut;
  datum_ = Data_.find( frame );
  if (datum_ == Data_.end()) return 1;
  //mprintf("DEBUG: Double frame %i is %lf\n",frame,(*it).second);
  *value = (*datum_).second;
  return 0;
}

// DataSet_double::Dval()
double DataSet_double::Dval(int idx) {
  double val;
  if (Get(&val,idx)) return 0;
  return val;
}

// DataSet_double::FrameIsEmpty()
int DataSet_double::FrameIsEmpty(int frame) {
  datum_ = Data_.find( frame );
  if (datum_ == Data_.end()) return 1;
  return 0;
}

// DataSet_double::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_double::WriteBuffer(CharBuffer &cbuffer, int frame) {
  double dval;
  datum_ = Data_.find( frame );
  if (datum_ == Data_.end())
    dval = 0.0;
  else
    dval = (*datum_).second;
  cbuffer.WriteDouble(data_format_, dval);
}

// DataSet_double::Width()
int DataSet_double::Width() {
  return (width_ + leadingSpace_);
}

// DataSet_double::Sync()
/** Since it seems to be very difficult (or impossible) to define Classes
  * as MPI datatypes, first non-master threads need to convert their maps
  * into 2 arrays, an int array containing frame #s and a double array
  * containing mapped values. These arrays are then sent to the master,
  * where they are converted pairs and inserted into the master map.
  */
int DataSet_double::Sync() {
  int rank, i, dataSize;
  int *Frames;
  double *Values;

  if (worldsize==1) return 0;

  for (rank = 1; rank < worldsize; rank ++) {
    // Get size of map on rank 
    if (worldrank>0) {
      // NOTE: current should be equal to size(). Check for now
      /*rprintf("DataSet_double syncing. current=%i, size=%u\n",
              current_, Data_.size());
      if (current_ != (int) Data_.size()) {
        rprintf("ERROR: current and map size are not equal.\n");
        return 1;
      }*/
      dataSize = (int)Data_.size();
    }

    // Send size of map on rank to master, allocate arrays on rank and master
    parallel_sendMaster(&dataSize, 1, rank, 0);
    rprintf("DataSet_double allocating %i for send/recv\n",dataSize);
    Frames = new int[ dataSize ];
    Values = new double[ dataSize ];
      
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
        Data_.insert( pair<int,double>( Frames[i], Values[i] ) );
    }

    // Free arrays
    delete[] Frames;
    delete[] Values;
  } // End loop over ranks>0

  return 0;
}

