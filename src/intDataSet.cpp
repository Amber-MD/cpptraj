// intDataSet
#include <cstdio>
#include <cstdlib>
#include "intDataSet.h"
#include "PtrajMpi.h"
#include "CpptrajStdio.h"
using namespace std;
/*
 * intDataSet::Add()
 * Insert data vIn at frame.
 */
void intDataSet::Add(int frame, void *vIn) {
  int *value;

  value = (int*) vIn;

  // Always insert at the end
  it=Data.end();
  Data.insert( it, pair<int,int>(frame, *value) );
  current++;
}

/*
 * intDataSet::isEmpty()
 */
int intDataSet::isEmpty(int frame) {
  it = Data.find( frame );
  if (it == Data.end()) return 1;
  return 0;
}

/*
 * intDataSet::Write()
 * Write data at frame to buffer. If no data for frame write 0.
 * Return position in buffer after write.
 */
char *intDataSet::Write(char *buffer, int frame) {

  it = Data.find( frame );
  if (it == Data.end()) 
    //sprintf(buffer," %12s","NoData");
    sprintf(buffer," %12i", 0);
  else 
    sprintf(buffer," %12i",(*it).second);
  return (buffer + 13);
}

/*
 * intDataSet::Width()
 */
int intDataSet::Width() {
  return 13;
}

/*
 * intDataSet::Sync()
 * Since it seems to be very difficult (or impossible) to define Classes
 * as MPI datatypes, first non-master threads need to convert their maps
 * into 2 arrays, an int array containing frame #s and another int array
 * containing mapped values. These arrays are then sent to the master,
 * where they are converted pairs and inserted into the master map.
 */
int intDataSet::Sync() {
  int rank, i, dataSize;
  int *Frames;
  int *Values;

  if (worldsize==1) return 0;

  for (rank = 1; rank < worldsize; rank ++) {
    // Get size of map on rank 
    if (worldrank>0) {
      // NOTE: current should be equal to size(). Check for now
      rprintf( "intDataSet syncing. current=%i, size=%u\n",
              current, Data.size());
      if (current != (int) Data.size()) {
        rprintf("ERROR: current and map size are not equal.\n");
        return 1;
      }
      dataSize = current;
    }

    // Send size of map on rank to master, allocate arrays on rank and master
    parallel_sendMaster(&dataSize, 1, rank, 0);
    rprintf("intDataSet allocating %i for send/recv\n",dataSize);
    Frames = (int*) malloc(dataSize * sizeof(int));
    Values = (int*) malloc(dataSize * sizeof(double));
      
    // On non-master convert map to 2 int arrays.
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
    parallel_sendMaster(Values, dataSize, rank, 0);

    // On master convert arrays to pairs and insert to master map
    if (worldrank==0) {
      for (i=0; i < dataSize; i++) 
        Data.insert( pair<int,int>( Frames[i], Values[i] ) );
    }

    // Free arrays
    free(Frames);
    free(Values);
  } // End loop over ranks>0

  return 0;
}

