#ifndef INC_DATASET_FLOAT_H
#define INC_DATASET_FLOAT_H
/// Class: DataSet_float 
/// Use the C++ STL map class instead of a straight array of floats. This
/// will allow Y values with non-consecutive X values to be stored. This is the
/// case e.g. when an action is not active for a certain part of the trajectory
/// read when it is not valid for the current parmtop.
#include <map>
#include "DataSet.h"
class DataSet_float : public DataSet {
    std::map<int,float> Data;
    std::map<int,float>::iterator datum;
  public:
    DataSet_float();

    int Xmax();
    int isEmpty(int);
    void Add( int, void * );
    int Get(void *, int);
    void WriteBuffer(CharBuffer&,int);
    int Width();
    int Sync();
};
#endif
