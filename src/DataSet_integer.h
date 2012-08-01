#ifndef INC_DATASET_INTEGER_H
#define INC_DATASET_INTEGER_H
#include <vector>
#include "DataSet.h"
// Class: DataSet_integer
/// Hold an array of integer values.
/** Actually 2 arrays; one for data and one for frame indices. This allows 
  * Y values with non-consecutive X values to be stored, which can happen 
  * e.g. when an action is not active for a certain trajectory because it 
  * is not valid for that topology.
  */
// TODO: Currently Resize and the [] operator are only used by Clustering
//       to create the cnumvtime array. Should these functions be
//       in all atomic datasets?
class DataSet_integer : public DataSet {
  public:
    DataSet_integer();

    void Resize(int);
    int& operator[](int idx) { return Data_[idx]; }

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
    typedef std::vector<int> DType;
    typedef std::vector<int> IType;
    DType Data_;
    DType::iterator datum_;
    IType Frames_;
    IType::iterator frame_;
};
#endif
