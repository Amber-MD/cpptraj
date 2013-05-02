#ifndef INC_DATASET_DOUBLE_H
#define INC_DATASET_DOUBLE_H
#include <vector>
#include "DataSet_1D.h"
// Class: DataSet_double
/// Hold an array of double values.
class DataSet_double : public DataSet_1D {
  public:
    DataSet_double() : DataSet_1D(DOUBLE, 12, 4) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_double();}
    double& operator[](size_t idx) { return Data_[idx];         }
    /// Make set size sizeIn, all values set to 0.0.
    void Resize(size_t sizeIn)     { Data_.resize(sizeIn, 0.0); }
    // ----- DataSet functions -------------------
    size_t Size()            const { return Data_.size();       }
    int Sync();
    void Info()              const { return;                    }
    // ----- DataSet_1D functions ----------------
    int Allocate1D(size_t);
    void Add( size_t, const void* );
    double Dval(size_t idx)  const { return Data_[idx];         }
    void WriteBuffer(CpptrajFile&, size_t) const;
    // -------------------------------------------
    typedef std::vector<double>::iterator iterator;
    iterator begin() { return Data_.begin(); }
    iterator end()   { return Data_.end();   }
  private:
    std::vector<double> Data_;
};
#endif
