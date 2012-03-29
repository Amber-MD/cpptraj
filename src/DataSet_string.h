#ifndef INC_DATASET_STRING_H
#define INC_DATASET_STRING_H
#include <map>
#include <string>
#include "DataSet.h"
// Class: DataSet_string
/// Hold an array of strings.
/** Use the C++ STL map and string classes instead of a straight array of 
  * char*. This will allow Y values with non-consecutive X values to be 
  * stored (e.g. when an action is not active for a certain part of the 
  * analysis) and puts the memory management burden on STL.
  */
class DataSet_string : public DataSet {
  public:
    DataSet_string();

    int Xmax();
    int Size();
    int FrameIsEmpty(int);
    void Add( int, void * );
    void WriteBuffer(CharBuffer&, int);
    int Width();
    int Sync();
  private:
    std::map<int,std::string> Data_;
    std::map<int,std::string>::iterator datum_;
};
#endif
