// mapDataSet
#include <cstdio>
#include <cstdlib>
#include "mapDataSet.h"
#include "PtrajMpi.h"
#include "CpptrajStdio.h"
using namespace std;
/*
 * mapDataSet::Add()
 * Insert data vIn at frame.
 */
void mapDataSet::Add(int frame, void *vIn) {
  double *value;

  value = (double*) vIn;

  // Always insert at the end
  it=Data.end();
  Data.insert( it, pair<int,double>(frame, *value) );
  current++;
}

/*
 * mapDataSet::isEmpty()
 */
int mapDataSet::isEmpty(int frame) {
  it = Data.find( frame );
  if (it == Data.end()) return 1;
  return 0;
}

/*
 * mapDataSet::Write()
 * Write data at frame to buffer, return 0. 
 * If no data for frame write 0.0, return 1.
 */
void mapDataSet::Write(char *buffer, int frame) {

  it = Data.find( frame );
  if (it == Data.end()) 
    //sprintf(buffer," %12s","NoData");
    sprintf(buffer," %12.4lf", 0.0);
  else 
    sprintf(buffer," %12.4lf",(*it).second);
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
  int rank, i, dataSize;
  int *Frames;
  double *Values;

  if (worldsize==1) return 0;

  for (rank = 1; rank < worldsize; rank ++) {
    // Get size of map on rank 
    if (worldrank>0) {
      // NOTE: current should be equal to size(). Check for now
      rprintf("mapDataSet syncing. current=%i, size=%u\n",
              current, Data.size());
      if (current != (int) Data.size()) {
        rprintf("ERROR: current and map size are not equal.\n");
        return 1;
      }
      dataSize = current;
    }

    // Send size of map on rank to master, allocate arrays on rank and master
    parallel_sendMaster(&dataSize, 1, rank, 0);
    rprintf("mapDataSet allocating %i for send/recv\n",dataSize);
    Frames = (int*) malloc(dataSize * sizeof(int));
    Values = (double*) malloc(dataSize * sizeof(double));
      
    // On non-master convert map to int and double arrays.
    if (worldrank > 0) {
      i=0;
      for ( it = Data.begin(); it != Data.end(); it++ ) {
        Frames[i]=(*it).first;
        Values[i]=(*it).second;
        i++;
      }
    }

    // Send arrays to master
    parallel_sendMaster(Frames, dataSize, rank, 0);
    parallel_sendMaster(Values, dataSize, rank, 1);

    // On master convert arrays to pairs and insert to master map
    if (worldrank==0) {
      for (i=0; i < dataSize; i++) 
        Data.insert( pair<int,double>( Frames[i], Values[i] ) );
    }

    // Free arrays
    free(Frames);
    free(Values);
  } // End loop over ranks>0

  return 0;
}

