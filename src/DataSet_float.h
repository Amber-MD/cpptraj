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
    /// Resize, filling extra values with given value
    void Resize(size_t sizeIn, float f)  { Data_.resize(sizeIn, f);   }
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
    size_t MemUsageInBytes() const { return Data_.size() * sizeof(float); }
    int MemAlloc(SizeArray const&);
    void CopyBlock(size_t, const DataSet*, size_t, size_t);
    // ----- DataSet_1D functions ----------------
    double Dval(size_t idx)         const { return (double)Data_[idx]; }
    double Xcrd(size_t idx)         const { return Dim(0).Coord(idx);  }
    const void* VoidPtr(size_t idx) const { return (void*)(&(Data_[0])+idx); }
    void SetY(size_t idx, double y)       { Data_[idx] = (float)y; }
    // -------------------------------------------
    float* Ptr()                         { return &(Data_[0]);        }
  private:
    std::vector<float> Data_;
};
#endif
