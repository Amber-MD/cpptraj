#ifndef INC_DATASET_STRING_H
#define INC_DATASET_STRING_H
#include <vector>
#include <string>
#include "DataSet.h"
/// Hold an array of string values.
class DataSet_string : public DataSet {
  public:
    DataSet_string() : DataSet(STRING, GENERIC, TextFormat(TextFormat::STRING, 1), 1) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_string();}
    std::vector<std::string> const& Data()          { return Data_;      }
    std::string& operator[](size_t idx)             { return Data_[idx]; }
    std::string const& operator[](size_t idx) const { return Data_[idx]; }

    void operator=(std::vector<std::string> const& rhs) { Data_ = rhs;}
    void AddElement(std::string const& s){ Data_.push_back( s );      }
    /// Make set size sizeIn, all values set to blank
    void Resize(size_t sizeIn)           { Data_.resize(sizeIn, "");  }
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
    size_t MemUsageInBytes() const;
  private:
    std::vector<std::string> Data_;
};
#endif
