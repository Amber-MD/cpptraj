#ifndef INC_DATASET_INTEGER_H
#define INC_DATASET_INTEGER_H
#include <vector>
#include "DataSet_1D.h"
// Class: DataSet_integer
/// Hold an array of integer values.
class DataSet_integer : public DataSet_1D {
  public:
    DataSet_integer() : DataSet_1D(INTEGER, 12, 0) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_integer();}
    int& operator[](size_t idx)    { return Data_[idx];         }
    /// Make set size sizeIn, all values set to 0.0.
    void Resize(size_t sizeIn)     { Data_.resize(sizeIn, 0);   }
    // ----- DataSet functions -------------------
    size_t Size()            const { return Data_.size();       }
    int Sync();
    void Info()              const { return;                    }
    // ----- DataSet_1D functions ----------------
    int Allocate1D(size_t);
    void Add( size_t, const void* );
    double Dval(size_t idx)  const { return (double)Data_[idx]; }
    void WriteBuffer(CpptrajFile&, size_t) const;
    // -------------------------------------------
    typedef std::vector<int>::iterator iterator;
    iterator begin() { return Data_.begin(); }
    iterator end()   { return Data_.end();   }
  private:
    std::vector<int> Data_;
};
#endif
