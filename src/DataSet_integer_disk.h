#ifndef INC_DATASET_INTEGER_DISK_H
#define INC_DATASET_INTEGER_DISK_H
#include "DataSet_integer.h"
class DataSet_integer_disk : public DataSet_integer {
  public:
    DataSet_integer_disk();
    static DataSet* Alloc() { return (DataSet*)new DataSet_integer_disk(); }
    // ----- DataSet functions -------------------
    size_t Size() const { return nvals_; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
#   endif
    void Info() const;
    int Allocate(SizeArray const&);
    void Add( size_t, const void* );
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    int Append(DataSet*);
    // ----- DataSet_1D functions ----------------
    double Dval(size_t) const;
    /// This function is invalid for DataSet_integer_disk
    const void* VoidPtr(size_t) const;
    // ----- DataSet_integer functions -----------
    int& operator[](size_t);
    int  operator[](size_t) const;
    void AddElement(int);
    void Resize(size_t);
    void AddVal(size_t, int);
    //typedef std::vector<int>::iterator iterator;
    //iterator begin()                  { return Data_.begin();      }
    //iterator end()                    { return Data_.end();        }
    //int* Ptr()                        { return &(Data_[0]);        }
  private:
    int ncid_;
    size_t start_[1]; ///< Hold current index
    size_t count_[1]; ///< Current size to read/write
    unsigned int nvals_; ///< Total number of values in data set
};
#endif
