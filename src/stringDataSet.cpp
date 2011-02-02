// stringDataSet
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "stringDataSet.h"
#include "PtrajMpi.h"
#include "CpptrajStdio.h"
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
 * Write data at frame to buffer. If no data for frame write NoData.
 * Return position in buffer after write.
 */
char *stringDataSet::Write(char *buffer, int frame) {

  it = Data.find( frame );
  if (it == Data.end()) { 
    sprintf(buffer," %s", "NoData");
    return (buffer + 7);
  } else 
    sprintf(buffer," %s",(*it).second.c_str());

  return (buffer + (*it).second.size() + 1);
}

/*
 * stringDataSet::Sync()
 * Since it seems to be very difficult (or impossible) to define Classes
 * as MPI datatypes, first non-master threads need to convert their maps
 * into 2 arrays, an int array containing frame #s and a char array
 * containing mapped values. An additional array of ints holding the size
 * of each string is also needed since the char array will be sent as
 * 1 big chunk. These arrays are then sent to the master, where they are 
 * converted to pairs and inserted into the master map.
 */
int stringDataSet::Sync() {
  int rank, i, dataSize, totalCharSize;
  int *Frames;
  char *Values;
  int *Sizes;
  char *ptr, oldChar;
  string str;

  if (worldsize==1) return 0;

  for (rank = 1; rank < worldsize; rank ++) {
    // Get size of map on rank 
    if (worldrank>0) {
      // NOTE: current should be equal to size(). Check for now
      rprintf( "stringDataSet syncing. current=%i, size=%u\n",
              current, Data.size());
      if (current != (int) Data.size()) {
        rprintf("ERROR: current and map size are not equal.\n");
        return 1;
      }
      dataSize = current;
    }

    // Send size of map on rank to master, allocate frame and sizes arrays on 
    // rank and master.
    parallel_sendMaster(&dataSize, 1, rank, 0);
    rprintf("stringDataSet allocating %i for send/recv\n",dataSize);
    Frames = (int*) malloc(dataSize * sizeof(int));
    Sizes = (int*) malloc(dataSize * sizeof(int));

    // On non-master get the size of each string in map and total size for
    // char array. Also get frame #s.
    if (worldrank>0) {
      i=0;
      totalCharSize=0;
      for (it = Data.begin(); it != Data.end(); it++) {
        Frames[i]=(*it).first;
        Sizes[i] = (*it).second.size();
        totalCharSize += (*it).second.size();
        i++;
      }
      // Add 1 to totalCharSize for NULL character
      totalCharSize++;
    }

    // Send total size of char array to master. allocate char array on rank and master
    parallel_sendMaster(&totalCharSize, 1, rank, 0);
    rprintf("stringDataSet allocating %i for char send/recv\n",totalCharSize);
    Values = (char*) malloc( totalCharSize * sizeof(char));
    strcpy(Values,"");

    // On non-master put all strings into giant char array
    if (worldrank > 0) {
      for (it = Data.begin(); it != Data.end(); it++) 
        strcat(Values, (*it).second.c_str());
    }

    // Send arrays to master
    parallel_sendMaster(Frames, dataSize, rank, 0);
    parallel_sendMaster(Sizes, dataSize, rank, 1);
    parallel_sendMaster(Values, totalCharSize, rank, 2);

    // On master convert arrays to pairs and insert to master map
    if (worldrank==0) {
      ptr = Values;
      for (i=0; i < dataSize; i++) { 
        // Recover string from giant char array
        oldChar = ptr[ Sizes[i] ];
        ptr[ Sizes[i] ] = '\0';
        str.assign(ptr);
        // Insert Frame/string pair
        Data.insert( pair<int,string>( Frames[i], str ) );
        // Advance pointer into giant char array
        ptr[ Sizes[i] ] = oldChar;
        ptr += Sizes[i];
      }
    }

    // Free arrays
    free(Frames);
    free(Sizes);
    free(Values);
  } // End loop over ranks>0

  return 0;
}

