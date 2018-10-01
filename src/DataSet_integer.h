#ifndef INC_DATASET_INTEGER_H
#define INC_DATASET_INTEGER_H
#include "DataSet_1D.h"
/// Base class for integer 1D data sets 
class DataSet_integer : public DataSet_1D {
  public:
    DataSet_integer() : DataSet_1D(INTEGER, TextFormat(TextFormat::INTEGER, 12)) {}
    virtual ~DataSet_integer() {} // Virtual bc inherited
    //static DataSet* Alloc() { return (DataSet*)new DataSet_integer();} TODO fix
    //virtual int& operator[](size_t) = 0;
    virtual void SetElement(size_t, int) = 0;
    virtual int  operator[](size_t) const = 0;
    virtual void AddElement(int) = 0;
    /// Make set size sizeIn, all values set to 0.0.
    virtual void Resize(size_t) = 0;
    /// Make set size sizeIn, all values set to val.
    virtual void Assign(size_t,int) = 0;
    virtual void AddVal(size_t, int) = 0;
    // ----- DataSet_1D functions ----------------
    double Xcrd(size_t idx)     const { return Dim(0).Coord(idx);  }
    // -------------------------------------------
    //typedef std::vector<int>::iterator iterator;
    //iterator begin()                  { return Data_.begin();      }
    //iterator end()                    { return Data_.end();        }
    //int* Ptr()                        { return &(Data_[0]);        }
};
#endif
