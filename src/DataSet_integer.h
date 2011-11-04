#ifndef INC_DATASET_INTEGER_H
#define INC_DATASET_INTEGER_H
/// Class: DataSet_integer
/// Use the C++ STL map class instead of a straight array of ints. This
/// will allow Y values with non-consecutive X values to be stored. This is the
/// case e.g. when an action is not active for a certain part of the analysis.
#include <map>
#include "DataSet.h"
class DataSet_integer : public DataSet {
    std::map<int,int> Data;
    std::map<int,int>::iterator it;
  public:
    DataSet_integer();

    int Xmax();
    int isEmpty(int);
    double Avg(double*);
    void Add( int, void * );
    int Get(void *, int);
    void WriteBuffer(CharBuffer&,int);
    int Width();
    int Sync();
};
#endif
