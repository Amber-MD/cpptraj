#ifndef INC_DATASET_INTEGER_MEM_H
#define INC_DATASET_INTEGER_MEM_H
#include <vector>
#include "DataSet_integer.h"
/// Hold an array of integer values in memory.
class DataSet_integer_mem : public DataSet_integer {
  public:
    DataSet_integer_mem() {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_integer_mem();}
    // ----- DataSet_integer functions -----------
    void SetElement(size_t idx, int val) { Data_[idx] = val; }
    int  operator[](size_t idx) const { return Data_[idx];         }
    void AddElement(int i)            { Data_.push_back( i );      }
    /// Make set size sizeIn, all values set to 0.0.
    void Resize(size_t sizeIn)        { Data_.resize(sizeIn, 0);   }
    /// Make set size sizeIn, all values set to val.
    void Assign(size_t sizeIn, int val) { Data_.resize(sizeIn, val); }
    inline void AddVal(size_t, int);
    // ----- DataSet_1D functions ----------------
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
    size_t MemUsageInBytes() const { return Data_.size() * sizeof(int); }
    // -------------------------------------------
    //typedef std::vector<int>::iterator iterator;
    //iterator begin()                  { return Data_.begin();      }
    //iterator end()                    { return Data_.end();        }
    //int* Ptr()                        { return &(Data_[0]);        }
  private:
    std::vector<int> Data_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
void DataSet_integer_mem::AddVal(size_t frame, int ival) {
  if (frame < Data_.size())
    Data_[frame] += ival;
  else {
    if (frame > Data_.size()) Data_.resize( frame, 0 );
    Data_.push_back( ival );
  }
}
#endif
