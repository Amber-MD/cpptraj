// DataSet_string
#include <cstring>
#include "DataSet_string.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"
using namespace std;

// CONSTRUCTOR
DataSet_string::DataSet_string() {
  width_ = 1;
  dType_ = STRING;
  SetDataSetFormat(false);
}

// DataSet_string::Xmax(()
/** Return the maximum X value added to this set. By convention this is 
  * always the last value added.
  */
int DataSet_string::Xmax() {
  // If no data has been added return 0
  //if (current_==0) return 0;
  if (Data_.empty()) return 0;
  datum_ = Data_.end();
  --datum_;
  return ( (*datum_).first );
}

// DataSet_string::Size()
int DataSet_string::Size() {
  return (int)Data_.size();
}

// DataSet_string::Add()
/** Insert data vIn at frame. If the size of the input string is greater
  * than the current width, increase the width.
  * String expects char*
  */
void DataSet_string::Add(int frame, void *vIn) {
  char *value;
  string Temp;
  int strsize;

  value = (char*) vIn;
  Temp.assign(value);
  strsize = (int) Temp.size();
  if (strsize > width_) width_ = strsize;
  // Always insert at the end
  //it=Data.end();
  //Data.insert( it, pair<int,string>(frame, Temp) );
  Data_[frame]=Temp;
  //++current_;
}

// DataSet_string::FrameIsEmpty()
int DataSet_string::FrameIsEmpty(int frame) {
  datum_ = Data_.find( frame );
  if (datum_ == Data_.end()) return 1;
  return 0;
}

// DataSet_string::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_string::WriteBuffer(CharBuffer &cbuffer, int frame) {
  datum_ = Data_.find( frame );
  if (datum_ == Data_.end())
    cbuffer.WriteString(data_format_,"NoData");
  else
    cbuffer.WriteString(data_format_, (*datum_).second.c_str());
}

// DataSet_string::Width()
/** Return the width in characters necessary to print data from this dataset.
  * Width is set whenever data is added and is the size of the largest stored
  * string.
  */
int DataSet_string::Width() {
  return (width_ + leadingSpace_);
}

// DataSet_string::Sync()
/** Since it seems to be very difficult (or impossible) to define Classes
  * as MPI datatypes, first non-master threads need to convert their maps
  * into 2 arrays, an int array containing frame #s and a char array
  * containing mapped values. An additional array of ints holding the size
  * of each string is also needed since the char array will be sent as
  * 1 big chunk. These arrays are then sent to the master, where they are 
  * converted to pairs and inserted into the master map.
  */
int DataSet_string::Sync() {
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
      //rprintf( "DataSet_string syncing. current=%i, size=%u\n",
      //        current_, Data_.size());
      /*if (current != (int) Data.size()) {
        rprintf("ERROR: current and map size are not equal.\n");
        return 1;
      }*/
      dataSize = (int)Data_.size();
    }

    // Send size of map on rank to master, allocate frame and sizes arrays on 
    // rank and master.
    parallel_sendMaster(&dataSize, 1, rank, 0);
    rprintf("DataSet_string allocating %i for send/recv\n",dataSize);
    Frames = new int[ dataSize ];
    Sizes = new int[ dataSize ];

    // On non-master get the size of each string in map and total size for
    // char array. Also get frame #s.
    if (worldrank>0) {
      i=0;
      totalCharSize=0;
      for (datum_ = Data_.begin(); datum_ != Data_.end(); datum_++) {
        Frames[i]=(*datum_).first;
        Sizes[i] = (*datum_).second.size();
        totalCharSize += (*datum_).second.size();
        i++;
      }
      // Add 1 to totalCharSize for NULL character
      totalCharSize++;
    }

    // Send total size of char array to master. allocate char array on rank and master
    parallel_sendMaster(&totalCharSize, 1, rank, 0);
    rprintf("DataSet_string allocating %i for char send/recv\n",totalCharSize);
    Values = new char[ totalCharSize ];
    strcpy(Values,"");

    // On non-master put all strings into giant char array
    if (worldrank > 0) {
      for (datum_ = Data_.begin(); datum_ != Data_.end(); datum_++) 
        strcat(Values, (*datum_).second.c_str());
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
        Data_.insert( pair<int,string>( Frames[i], str ) );
        // Advance pointer into giant char array
        ptr[ Sizes[i] ] = oldChar;
        ptr += Sizes[i];
      }
    }

    // Free arrays
    delete[] Frames;
    delete[] Sizes;
    delete[] Values;
  } // End loop over ranks>0

  return 0;
}

