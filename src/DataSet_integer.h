#ifndef INC_DATASET_INTEGER_H
#define INC_DATASET_INTEGER_H
#include <vector>
#include "DataSet_1D.h"
// Class: DataSet_integer
/// Hold an array of integer values.
class DataSet_integer : public DataSet_1D {
  public:
    DataSet_integer() : DataSet_1D(INTEGER, TextFormat(TextFormat::INTEGER, 12)) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_integer();}
    int& operator[](size_t idx)       { return Data_[idx];         }
    int  operator[](size_t idx) const { return Data_[idx];         }
    void AddElement(int i)            { Data_.push_back( i );      }
    /// Make set size sizeIn, all values set to 0.0.
    void Resize(size_t sizeIn)        { Data_.resize(sizeIn, 0);   }
    inline void AddVal(size_t, int);
    // ----- DataSet functions -------------------
    size_t Size()               const { return Data_.size();       }
    int Sync(size_t, std::vector<int> const&);
    void Info()                 const { return;                    }
    int Allocate(SizeArray const&);
    void Add( size_t, const void* );
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    int Append(DataSet*);
    // ----- DataSet_1D functions ----------------
    double Dval(size_t idx)     const { return (double)Data_[idx]; }
    double Xcrd(size_t idx)     const { return Dim(0).Coord(idx);  } 
    // -------------------------------------------
    typedef std::vector<int>::iterator iterator;
    iterator begin()                  { return Data_.begin();      }
    iterator end()                    { return Data_.end();        }
    int* Ptr()                        { return &(Data_[0]);        }
  private:
    std::vector<int> Data_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
void DataSet_integer::AddVal(size_t frame, int ival) {
  if (frame < Data_.size())
    Data_[frame] += ival;
  else {
    if (frame > Data_.size()) Data_.resize( frame, 0 );
    Data_.push_back( ival );
  }
}
#endif
