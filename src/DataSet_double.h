#ifndef INC_DATASET_DOUBLE_H
#define INC_DATASET_DOUBLE_H
/// Class: DataSet_double
/// Use the C++ STL map class instead of a straight array of doubles. This
/// will allow Y values with non-consecutive X values to be stored. This is the
/// case e.g. when an action is not active for a certain part of the trajectory
/// read when it is not valid for the current parmtop.
#include <map>
#include "DataSet.h"
class DataSet_double : public DataSet {
    std::map<int,double> Data;
    std::map<int,double>::iterator it;
  public:
    DataSet_double();

    int Xmax();
    double Avg(double*);
    double Min();
    double Max();
    int isEmpty(int);
    void Add( int, void * );
    int Get(void *, int);
    void WriteBuffer(CharBuffer&, int);
    int Width();
    int Sync();
};
#endif
