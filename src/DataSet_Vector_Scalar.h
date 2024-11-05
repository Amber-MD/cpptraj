#ifndef INC_DATASET_VECTOR_SCALAR_H
#define INC_DATASET_VECTOR_SCALAR_H
#include "DataSet.h"
#include "Vec3.h"
#include <cstddef> // size_t
/// Hold XYZ value and associated scalar. 
class DataSet_Vector_Scalar : public DataSet {
  public:
    typedef std::vector<double> Darray;

    DataSet_Vector_Scalar();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Vector_Scalar(); }
    // ----- DataSet functions -------------------
    size_t Size()                                    const { return vecs_.size(); }
    void Info()                                      const { return; }
    int Allocate(SizeArray const&);
    inline void Add(size_t, const void*);
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    int Append(DataSet*);
    size_t MemUsageInBytes() const;
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
    /// Ensure each process has all frames
    int Bcast(Parallel::Comm const&);
#   endif
    // -------------------------------------------
    Vec3 const& Vec(unsigned int i) const { return vecs_[i]; }
    double Val(unsigned int i)      const { return vals_[i]; }
    const double* ValPtr()          const { return &vals_[0]; }

    void Resize(size_t s)           { vecs_.resize( s ); vals_.resize( s ); }
  
    Vec3& ModifyVec(unsigned int i) { return vecs_[i]; }
    double* ValPtr()                { return &vals_[0]; } 
  private:
    typedef std::vector<Vec3> Varray;

    Varray vecs_; ///< Hold XYZ values.
    Darray vals_; ///< Hold scalar values.
};

/** Add XYZ scalar to the DataSet */
void DataSet_Vector_Scalar::Add(size_t idx, const void* ptrIn) {
  const double* ptr = (const double*)ptrIn;
  vecs_.push_back( Vec3(ptr[0], ptr[1], ptr[2]) );
  vals_.push_back( ptr[3] );
}
#endif
