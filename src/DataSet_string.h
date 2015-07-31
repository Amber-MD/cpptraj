#ifndef INC_DATASET_STRING_H
#define INC_DATASET_STRING_H
#include <vector>
#include <string>
#include "DataSet.h"
/// Hold an array of string values.
class DataSet_string : public DataSet {
  public:
    DataSet_string() : DataSet(STRING, GENERIC, 1, 0, 1) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_string();}
    std::string& operator[](size_t idx)  { return Data_[idx];         }
    void operator=(std::vector<std::string> const& rhs) { Data_ = rhs;}
    void AddElement(std::string const& s){ Data_.push_back( s );      }
    void Append(std::vector<std::string> const&);
    /// Make set size sizeIn, all values set to blank
    void Resize(size_t sizeIn)           { Data_.resize(sizeIn, "");  }
    // ----- DataSet functions -------------------
    size_t Size()                  const { return Data_.size();       }
    int Sync();
    void Info()                    const { return;                    }
    int Allocate(SizeArray const&);
    void Add( size_t, const void* );
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
  private:
    std::vector<std::string> Data_;
};
#endif
