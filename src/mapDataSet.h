#ifndef INC_MAPDATASET_H
#define INC_MAPDATASET_H
/* Class: mapDataSet
 * Test of the C++ STL map class instead of a straight array of doubles. This
 * will allow Y values with non-consecutive X values to be stored. This is the
 * case e.g. when an action is not active for a certain part of the analysis.
 */
#include <map>
#include "DataSet.h"
//using namespace std;
class mapDataSet : public DataSet {
    std::map<int,double> Data;
    std::map<int,double>::iterator it;
  public:
    int isEmpty(int);
    void Add( int, void * );
    char *Write(char *, int);
    int Width();
    int Sync();
};
#endif
