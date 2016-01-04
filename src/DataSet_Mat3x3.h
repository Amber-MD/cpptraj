#ifndef INC_DATASET_MAT3X3_H
#define INC_DATASET_MAT3X3_H
#include "DataSet.h"
#include "Matrix_3x3.h"
/// Time series of Matrix_3x3 class.
class DataSet_Mat3x3 : public DataSet {
    typedef std::vector<Matrix_3x3> Marray;
  public:
    DataSet_Mat3x3() : DataSet(MAT3X3,GENERIC,TextFormat(TextFormat::DOUBLE,12,9,9),1) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_Mat3x3(); }
    // ----- DataSet functions -------------------
    size_t Size()                       const { return data_.size(); }
    int Sync(size_t, std::vector<int> const&);
    void Info()                         const { return;              }
    int Allocate(SizeArray const&);
    inline void Add(size_t, const void*);
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    int Append(DataSet*);
    // -------------------------------------------
    void AddMat3x3( Matrix_3x3 const& m) { data_.push_back( m ); }
    typedef Marray::const_iterator const_iterator;
    const_iterator begin() const { return data_.begin(); }
    const_iterator end()   const { return data_.end();   }
    typedef Marray::iterator iterator;
    iterator begin() { return data_.begin(); }
    iterator end()   { return data_.end();   }
    Matrix_3x3 const& operator[](int i) { return data_[i]; }
  private:
    Marray data_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
void DataSet_Mat3x3::Add(size_t frame, const void* vIn) {
  if (frame > data_.size())
    data_.resize(frame, Matrix_3x3(0.0));
  data_.push_back( Matrix_3x3( (const double*)vIn ) );
}
#endif
