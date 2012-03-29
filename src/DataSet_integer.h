#ifndef INC_DATASET_INTEGER_H
#define INC_DATASET_INTEGER_H
#include <map>
#include "DataSet.h"
// Class: DataSet_integer
/// Hold an array of integers.
/** Use the C++ STL map class instead of a straight array of ints. This
  * will allow Y values with non-consecutive X values to be stored. This is the
  * case e.g. when an action is not active for a certain part of the analysis.
  */
class DataSet_integer : public DataSet {
  public:
    DataSet_integer();

    void Begin();
    bool NextValue();
    double CurrentValue();

    int Xmax();
    int Size();
    int FrameIsEmpty(int);
    void Add( int, void * );
    int Get(void *, int);
    double Dval(int);
    void WriteBuffer(CharBuffer&,int);
    int Width();
    int Sync();
  private:
    std::map<int,int> Data_;
    std::map<int,int>::iterator datum_;
};
#endif
