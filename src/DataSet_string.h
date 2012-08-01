#ifndef INC_DATASET_STRING_H
#define INC_DATASET_STRING_H
#include <vector>
#include <string>
#include "DataSet.h"
// Class: DataSet_string
/// Hold an array of strings.
/** Actually 2 arrays; one for data and one for frame indices. This allows 
  * Y values with non-consecutive X values to be stored, which can happen 
  * e.g. when an action is not active for a certain trajectory because it 
  * is not valid for that topology.
  */
class DataSet_string : public DataSet {
  public:
    DataSet_string();

    int Allocate(int);
    int Xmax();
    int Size();
    int FrameIsEmpty(int);
    void Add( int, void * );
    void WriteBuffer(CpptrajFile&, int);
    int Sync();
  private:
    typedef std::vector<std::string> DType;
    typedef std::vector<int> IType;
    DType Data_;
    DType::iterator datum_;
    IType Frames_;
    IType::iterator frame_;
};
#endif
