#ifndef INC_DATASET_STRINGVAR_H
#define INC_DATASET_STRINGVAR_H
#include "DataSet.h"
#include <string>
/// Hold a single string variable 
class DataSet_StringVar : public DataSet {
  public:
    DataSet_StringVar();
    static DataSet* Alloc() { return (DataSet*)new DataSet_StringVar(); }
    // ----- DataSet functions -------------------
    size_t Size()                                    const { return 1; }
    void Info()                                      const { return; }
    int Allocate(SizeArray const&)                         { return 0; }
    void Add(size_t, const void*);
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    int Append(DataSet*)                                   { return 1; }
    size_t MemUsageInBytes() const;
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    // -------------------------------------------

    void append(std::string const& s) { var_.append(s); }
    void assign(std::string const& s) { var_.assign(s); }

    std::string const& Value() const { return var_; }
  private:
    std::string var_; ///< String variable
};
#endif
