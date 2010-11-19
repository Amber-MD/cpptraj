#ifndef INC_STRINGDATASET_H
#define INC_STRINGDATASET_H
/* Class: stringDataSet
 * Test of the C++ STL map class instead of a straight array of strings. This
 * will allow Y values with non-consecutive X values to be stored. This is the
 * case e.g. when an action is not active for a certain part of the analysis.
 */
#include <map>
#include <string>
#include "DataSet.h"
class stringDataSet : public DataSet {
    std::map<int,std::string> Data;
    std::map<int,std::string>::iterator it;
  public:
    int isEmpty(int);
    void Add( int, void *vIn );
    void Write(char *, int);
    //int Sync();
};
#endif
