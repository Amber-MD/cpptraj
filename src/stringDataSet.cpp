// stringDataSet
//#include <cstdlib>
#include <cstdio>
#include "stringDataSet.h"
//#include "PtrajMpi.h"
using namespace std;
/*
 * stringDataSet::Add()
 * Insert data vIn at frame.
 * String expects char*
 */
void stringDataSet::Add(int frame, void *vIn) {
  char *value;
  string Temp;

  value = (char*) vIn;
  Temp.assign(value);
  // Always insert at the end
  it=Data.end();
  Data.insert( it, pair<int,string>(frame, Temp) );
  current++;
}

/*
 * stringDataSet::isEmpty()
 */
int stringDataSet::isEmpty(int frame) {
  it = Data.find( frame );
  if (it == Data.end()) return 1;
  return 0;
}

/*
 * stringDataSet::Write()
 * Write data at frame to buffer and return 0. 
 * If no data for frame write NoData and return 1.
 */
void stringDataSet::Write(char *buffer, int frame) {

  it = Data.find( frame );
  if (it == Data.end()) 
    sprintf(buffer," %s", "NoData");
  else 
    sprintf(buffer," %s",(*it).second.c_str());
}

/*
 * stringDataSet::Sync()
 * Since it seems to be very difficult (or impossible) to define Classes
 * as MPI datatypes, first non-master threads need to convert their maps
 * into 2 arrays, an int array containing frame #s and a double array
 * containing mapped values. These arrays are then sent to the master,
 * where they are converted pairs and inserted into the master map.
 */
/*int stringDataSet::Sync() {
  int rank, i, dataSize;
  int *Frames;
  double *Values;

  if (worldsize==1) return 0;

  for (rank = 1; rank < worldsize; rank ++) {
    // Get size of map on rank 
    if (worldrank>0) {
      // NOTE: current should be equal to size(). Check for now
      rprintf(stdout, "mapDataSet syncing. current=%i, size=%u\n",
              current, Data.size());
      if (current != (int) Data.size()) {
        rprintf(stdout,"ERROR: current and map size are not equal.\n");
        return 1;
      }
      dataSize = current;
    }

    // Send size of map on rank to master, allocate arrays on rank and master
    parallel_sendMaster(&dataSize, 1, rank, 0);
    rprintf(stdout,"mapDataSet allocating %i for send/recv\n",dataSize);
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
*/
