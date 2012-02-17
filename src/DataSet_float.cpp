// DataSet_float
#include "DataSet_float.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"
using namespace std;

// CONSTRUCTOR
DataSet_float::DataSet_float() {
  width = 8;
  precision = 3;
  dType=FLOAT;
  SetDataSetFormat(false);
}

// DataSet_float::Begin()
void DataSet_float::Begin() {
  datum = Data.begin();
}

// DataSet_float::NextValue();
bool DataSet_float::NextValue() {
  datum++;
  if (datum == Data.end()) return false;
  return true;
}

// DataSet_float::CurrentValue()
double DataSet_float::CurrentValue() {
  return (double)(*datum).second;
}


/* DataSet_float::Xmax(()
 * Return the maximum X value added to this set. By convention this is 
 * always the last value added.
 */
int DataSet_float::Xmax() {
  // If no data has been added return 0
  if (current==0) return 0;
  datum=Data.end();
  datum--;
  return ( (*datum).first );
} 

/* DataSet_float::Add()
 * Insert data vIn at frame.
 */
void DataSet_float::Add(int frame, void *vIn) {
  float *value = (float*) vIn;

  // Always insert at the end
  //datum=Data.end();
  //Data.insert( datum, pair<int,float>(frame, *value) );
  Data[frame] = (*value);
  ++current;
}

/* DataSet_float::Get()
 * Get data at frame, put into vOut. Return 1 if no data at frame.
 */
int DataSet_float::Get(void *vOut, int frame) {
  float *value;
  
  if (vOut==NULL) return 1;
  //mprintf("DEBUG: Attempting to get float frame %i\n",frame);
  value = (float*) vOut;
  datum=Data.find( frame );
  if (datum == Data.end()) return 1;
  //mprintf("DEBUG: Float frame %i is %lf\n",frame,(*datum).second);
  *value = (*datum).second;
  return 0;
}

// DataSet_float::Dval()
double DataSet_float::Dval(int idx) {
  float val;
  if (Get(&val,idx)) return 0;
  return (double)val;
}

/* DataSet_float::isEmpty()
 */
int DataSet_float::isEmpty(int frame) {
  datum = Data.find( frame );
  if (datum == Data.end()) return 1;
  return 0;
}

/* DataSet_float::WriteBuffer()
 * Write data at frame to CharBuffer. If no data for frame write 0.0.
 */
void DataSet_float::WriteBuffer(CharBuffer &cbuffer, int frame) {
  double dval;
  datum = Data.find( frame );
  if (datum == Data.end())
    dval = 0.0;
  else
    dval = (double)(*datum).second;
  cbuffer.WriteDouble(data_format, dval);
}

/* DataSet_float::Width()
 */
int DataSet_float::Width() {
  return (width + leadingSpace);
}

/* DataSet_float::Sync()
 * Since it seems to be very difficult (or impossible) to define Classes
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
      rprintf("DataSet_float syncing. current=%i, size=%u\n",
              current, Data.size());
      if (current != (int) Data.size()) {
        rprintf("ERROR: current and map size are not equal.\n");
        return 1;
      }
      dataSize = current;
    }

    // Send size of map on rank to master, allocate arrays on rank and master
    parallel_sendMaster(&dataSize, 1, rank, 0);
    rprintf("DataSet_float allocating %i for send/recv\n",dataSize);
    Frames = new int[ dataSize ];
    Values = new float[ dataSize ];
      
    // On non-master convert map to int and double arrays.
    if (worldrank > 0) {
      i=0;
      for ( datum = Data.begin(); datum != Data.end(); datum++ ) {
        Frames[i]=(*datum).first;
        Values[i]=(*datum).second;
        i++;
      }
    }

    // Send arrays to master
    parallel_sendMaster(Frames, dataSize, rank, 0);
    parallel_sendMaster(Values, dataSize, rank, 1);

    // On master convert arrays to pairs and insert to master map
    if (worldrank==0) {
      for (i=0; i < dataSize; i++) 
        Data.insert( pair<int,float>( Frames[i], Values[i] ) );
    }

    // Free arrays
    delete[] Frames;
    delete[] Values;
  } // End loop over ranks>0

  return 0;
}

