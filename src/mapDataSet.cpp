// mapDataSet
#include <cstdio>
#include <cstdlib>
#include "mapDataSet.h"
#include "PtrajMpi.h"
#include "CpptrajStdio.h"
using namespace std;

// CONSTRUCTOR
mapDataSet::mapDataSet() {
  width = 12;
  precision = 4;
  totalwidth = 39;
  dType=MAP;
  setFormatString();
}

/* mapDataSet::Xmax(()
 * Return the maximum X value added to this set. By convention this is 
 * always the last value added.
 */
int mapDataSet::Xmax() {
  // If no data has been added return 0
  if (current==0) return 0;
  return xData.back();
} 

/* mapDataSet::Add()
 * Insert data vIn. Expect an array of size 3. frame is not used in 
 * this case.
 */
void mapDataSet::Add(int frame, void *vIn) {
  double *value;

  value = (double*) vIn;

  // Always insert at the end
  xData.push_back(value[0]);
  yData.push_back(value[1]);
  zData.push_back(value[2]);

  current++;
}

/* mapDataSet::Get()
 * Get data at frame, put into vOut. Return 1 if no data at frame.
 */
int mapDataSet::Get(void *vOut, int frame) {
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

/* mapDataSet::isEmpty()
 * By definition no part of the map can be empty unless frame is out of bounds.
 */
int mapDataSet::isEmpty(int frame) {
  if (frame < 0 || frame >= (int)xData.size()) return 1;
  //it = Data.find( frame );
  //if (it == Data.end()) return 1;
  return 0;
}

/* mapDataSet::Write()
 * Write data at frame to buffer. If no data for frame write 0.0.
 * Return position in buffer after write.
 */
char *mapDataSet::Write(char *buffer, int frame) {

  if (isEmpty(frame)) 
    //sprintf(buffer," %12s","NoData");
    sprintf(buffer, format, 0.0, 0.0, 0.0);
  else 
    sprintf(buffer, format,xData[frame],yData[frame],zData[frame]);
  return (buffer + totalwidth);
}

/* mapDataSet::Width()
 */
int mapDataSet::Width() {
  return (totalwidth);
}

/*
 * mapDataSet::Sync()
 * Since it seems to be very difficult (or impossible) to define Classes
 * as MPI datatypes, first non-master threads need to convert their maps
 * into 2 arrays, an int array containing frame #s and a double array
 * containing mapped values. These arrays are then sent to the master,
 * where they are converted pairs and inserted into the master map.
 */
int mapDataSet::Sync() {
  int rank, idx, dataSize;
  double *Values;

  if (worldsize==1) return 0;

  for (rank = 1; rank < worldsize; rank ++) {
    // Get size of map on rank 
    if (worldrank>0) {
      // NOTE: current should be equal to size(). Check for now
      // Should really check all three datasets
      rprintf("mapDataSet syncing. current=%i, size=%u\n",
              current, xData.size());
      if (current != (int) xData.size()) {
        rprintf("ERROR: current and map size are not equal.\n");
        return 1;
      }
      dataSize = current;
    }

    // Send size of map on rank to master, allocate arrays on rank and master
    parallel_sendMaster(&dataSize, 1, rank, 0);
    rprintf("mapDataSet allocating %i for send/recv\n",dataSize);
    Values = (double*) malloc(3 * dataSize * sizeof(double));
      
    // On non-master convert map arrays to giant double array.
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
    free(Values);
  } // End loop over ranks>0

  return 0;
}

