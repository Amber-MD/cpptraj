// DataSet_XYZ
#include "DataSet_XYZ.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"
using namespace std;

// CONSTRUCTOR
DataSet_XYZ::DataSet_XYZ() {
  width = 12;
  precision = 4;
  totalwidth = 38;
  dType=XYZ;
  SetDataSetFormat(false);
}

// DataSet_XYZ::Xmax(()
/** Return the maximum X value added to this set. By convention this is 
  * always the last value added.
  */
int DataSet_XYZ::Xmax() {
  // If no data has been added return 0
  if (current==0) return 0;
  return xData.back();
} 

// DataSet_XYZ::Add()
/** Insert data vIn. Expect an array of size 3. frame is not used in 
  * this case.
  */
void DataSet_XYZ::Add(int frame, void *vIn) {
  double *value;

  value = (double*) vIn;

  // Always insert at the end
  xData.push_back(value[0]);
  yData.push_back(value[1]);
  zData.push_back(value[2]);

  current++;
}

// DataSet_XYZ::Get()
/** Get data at frame, put into vOut. 
  * \return 1 if no data at frame.
  */
int DataSet_XYZ::Get(void *vOut, int frame) {
  double *value;
  
  if (vOut==NULL) return 1;
  //mprintf("DEBUG: Attempting to get double frame %i\n",frame);
  value = (double*) vOut;
  if (frame < 0 || frame >= (int)xData.size()) return 1;
  value[0] = xData[frame];
  value[1] = yData[frame];
  value[2] = zData[frame];
  //mprintf("DEBUG: Double frame %i is %lf %lf %lf\n",frame,value[0],value[1],value[2]);
  return 0;
}

// DataSet_XYZ::isEmpty()
/** By definition no part of the map can be empty unless frame is out of bounds.
  */
int DataSet_XYZ::isEmpty(int frame) {
  if (frame < 0 || frame >= (int)xData.size()) return 1;
  //it = Data.find( frame );
  //if (it == Data.end()) return 1;
  return 0;
}

// DataSet_XYZ::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_XYZ::WriteBuffer(CharBuffer &cbuffer, int frame) {
  double darray[3];
  if (isEmpty(frame)) {
    darray[0] = 0;
    darray[1] = 0;
    darray[2] = 0;
  } else {
    darray[0] = xData[frame];
    darray[1] = yData[frame];
    darray[2] = zData[frame];
  }
  cbuffer.WriteDoubleXYZ(data_format, darray);
}


// DataSet_XYZ::Width()
int DataSet_XYZ::Width() {
  return (totalwidth + leadingSpace);
}

// DataSet_XYZ::Sync()
/** Since it seems to be very difficult (or impossible) to define Classes
  * as MPI datatypes, first non-master threads combine their arrays into
  * 1 giant array, which is then sent to master. 
  * NOTE: Should just be 3 separate sends? 
  */
int DataSet_XYZ::Sync() {
  int rank, idx, dataSize;
  double *Values;

  if (worldsize==1) return 0;

  for (rank = 1; rank < worldsize; rank ++) {
    // Get size of map on rank 
    if (worldrank>0) {
      // NOTE: current should be equal to size(). Check for now
      // Should really check all three datasets
      rprintf("DataSet_XYZ syncing. current=%i, size=%u\n",
              current, xData.size());
      if (current != (int) xData.size()) {
        rprintf("ERROR: current and map size are not equal.\n");
        return 1;
      }
      dataSize = current;
    }

    // Send size of map on rank to master, allocate arrays on rank and master
    parallel_sendMaster(&dataSize, 1, rank, 0);
    rprintf("DataSet_XYZ allocating %i for send/recv\n",dataSize);
    Values = new double[ 3 * dataSize ];
      
    // On non-master convert map arrays to giant double array.
    // NOTE: Could use memcopys eventually
    if (worldrank > 0) {
      idx=0;
      for ( int N = 0; N < dataSize; N++)
        Values[idx++]=xData[N];
      for ( int N = 0; N < dataSize; N++)
        Values[idx++]=yData[N];
      for ( int N = 0; N < dataSize; N++)
        Values[idx++]=zData[N];
    }

    // Send giant array to master
    parallel_sendMaster(Values, 3 * dataSize, rank, 1);

    // On master convert giant arrays to three arrays and insert to master map
    if (worldrank==0) {
      idx=0;
      for (int N=0; N < dataSize; N++)
        xData.push_back(Values[idx++]);
      for (int N=0; N < dataSize; N++)
        yData.push_back(Values[idx++]);
      for (int N=0; N < dataSize; N++)
        zData.push_back(Values[idx++]);
    }

    // Free arrays
    delete[] Values;
  } // End loop over ranks>0

  return 0;
}

