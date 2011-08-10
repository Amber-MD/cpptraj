#ifndef INC_DOUBLEDATASET_H
#define INC_DOUBLEDATASET_H
/// Class: doubleDataSet
/// Use the C++ STL map class instead of a straight array of doubles. This
/// will allow Y values with non-consecutive X values to be stored. This is the
/// case e.g. when an action is not active for a certain part of the trajectory
/// read when it is not valid for the current parmtop.
#include <map>
#include "DataSet.h"
class doubleDataSet : public DataSet {
    std::map<int,double> Data;
    std::map<int,double>::iterator it;
  public:
    doubleDataSet();

    int Xmax();
    int isEmpty(int);
    void Add( int, void * );
    int Get(void *, int);
    char *Write(char *, int);
    int Width();
    int Sync();
};
#endif
