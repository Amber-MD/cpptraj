#ifndef INC_DATASET_FLOAT_H
#define INC_DATASET_FLOAT_H
#include <vector>
#include "DataSet_1D.h"
/// Hold an array of float values.
class DataSet_float : public DataSet_1D {
  public:
    DataSet_float() : DataSet_1D(FLOAT, TextFormat(TextFormat::DOUBLE, 8, 3)) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_float();}
    float& operator[](size_t idx)        { return Data_[idx];         }
    float  operator[](size_t idx)  const { return Data_[idx];         }
    void AddElement(float f)             { Data_.push_back( f );      }
    /// Make set size sizeIn, all values set to 0.0.
    void Resize(size_t sizeIn)           { Data_.resize(sizeIn, 0.0); }
    // ----- DataSet functions -------------------
    size_t Size()                  const { return Data_.size();       }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
#   endif
    void Info()                    const { return;                    }
    int Allocate(SizeArray const&);
    void Add( size_t, const void* );
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    int Append(DataSet*);
    // ----- DataSet_1D functions ----------------
    double Dval(size_t idx)        const { return (double)Data_[idx]; }
    double Xcrd(size_t idx)        const { return Dim(0).Coord(idx);  }
  private:
    std::vector<float> Data_;
};
#endif
