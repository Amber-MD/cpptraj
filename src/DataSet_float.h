#ifndef INC_DATASET_FLOAT_H
#define INC_DATASET_FLOAT_H
#include <map>
#include "DataSet.h"
// Class: DataSet_float
/// Hold an array of floats. 
/** Use the C++ STL map class instead of a straight array of floats. This
  * will allow Y values with non-consecutive X values to be stored. This is the
  * case e.g. when an action is not active for a certain part of the trajectory
  * read when it is not valid for the current parmtop.
  */
class DataSet_float : public DataSet {
  public:
    DataSet_float();

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
    std::map<int,float> Data_;
    std::map<int,float>::iterator datum_;
};
#endif
