#ifndef INC_DATASET_UNSIGNEDINT_MEM_H
#define INC_DATASET_UNSIGNEDINT_MEM_H
#include <vector>
#include "DataSet_1D.h"
/// Hold an array of unsigned integer values in memory.
class DataSet_unsignedInt : public DataSet_1D {
  public:
    DataSet_unsignedInt() : DataSet_1D(UNSIGNED_INTEGER, TextFormat(TextFormat::UNSIGNED, 12)) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_unsignedInt();}
    // ----- DataSet_integer functions -----------
    void SetElement(size_t idx, unsigned int val) { Data_[idx] = val; }
    unsigned int operator[](size_t idx) const { return Data_[idx];         }
    void AddElement(unsigned int i)            { Data_.push_back( i );      }
    /// Make set size sizeIn, all values set to 0.0.
    void Resize(size_t sizeIn)        { Data_.resize(sizeIn, 0);   }
    /// Make set size sizeIn, all values set to val.
    void Assign(size_t sizeIn, unsigned int val) { Data_.resize(sizeIn, val); }
    inline void AddVal(size_t, unsigned int);
    // ----- DataSet_1D functions ----------------
    double Xcrd(size_t idx)     const { return Dim(0).Coord(idx);  }
    void SetY(size_t idx, double y) { SetElement(idx, (unsigned int)y); }
    double Dval(size_t idx)     const { return (double)Data_[idx]; }
    const void* VoidPtr(size_t idx) const { return (void*)(&(Data_[0])+idx); }
    // ----- DataSet functions -------------------
    size_t Size()               const { return Data_.size();       }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
    int Recv(size_t, unsigned int, int, int, int, Parallel::Comm const&);
    int Send(int, int, Parallel::Comm const&) const;
#   endif
    void Info()                 const { return;                    }
    int Allocate(SizeArray const&);
    void Add( size_t, const void* );
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    int Append(DataSet*);
    size_t MemUsageInBytes() const { return Data_.size() * sizeof(unsigned int); }
    int MemAlloc(SizeArray const&);
    void CopyBlock(size_t, const DataSet*, size_t, size_t);
    // -------------------------------------------
    //typedef std::vector<int>::iterator iterator;
    //iterator begin()                  { return Data_.begin();      }
    //iterator end()                    { return Data_.end();        }
    //int* Ptr()                        { return &(Data_[0]);        }
  private:
    std::vector<unsigned int> Data_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
void DataSet_unsignedInt::AddVal(size_t frame, unsigned int ival) {
  if (frame < Data_.size())
    Data_[frame] += ival;
  else {
    if (frame > Data_.size()) Data_.resize( frame, 0 );
    Data_.push_back( ival );
  }
}
#endif
