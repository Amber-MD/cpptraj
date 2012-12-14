#ifndef INC_DATASET_DOUBLE_H
#define INC_DATASET_DOUBLE_H
#include <vector>
#include "DataSet.h"
// Class: DataSet_double
/// Hold an array of double values.
/** Actually 2 arrays; one for data and one for frame indices. This allows 
  * Y values with non-consecutive X values to be stored, which can happen 
  * e.g. when an action is not active for a certain trajectory because it 
  * is not valid for that topology.
  */
class DataSet_double : public DataSet {
  public:
    DataSet_double();
    double& operator[](int idx) { return Data_[idx]; }
    void Resize(int);

    int Allocate(int);
    int Xmax();
    int Size();
    int FrameIsEmpty(int);
    void Add( int, void * );
    double CurrentDval();
    double Dval(int);
    void WriteBuffer(CpptrajFile&, int);
    int Sync();
  private:
    typedef std::vector<double> DType;
    typedef std::vector<int> IType;
    DType Data_;
    DType::iterator datum_;
    IType Frames_;
    IType::iterator frame_;
};
#endif
