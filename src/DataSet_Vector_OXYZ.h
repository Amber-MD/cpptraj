#ifndef INC_DATASET_VECTOR_OXYZ_H
#define INC_DATASET_VECTOR_OXYZ_H
#include "DataSet_Vector.h"
/// Vector data set with origins
class DataSet_Vector_OXYZ : public DataSet_Vector {
  public:
    DataSet_Vector_OXYZ();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Vector_OXYZ();}
    // ----- DataSet functions -------------------
    size_t MemUsageInBytes() const;
    int Allocate(SizeArray const&);
    void WriteBuffer(CpptrajFile&,SizeArray const&) const;
    inline void Add(size_t, const void*);
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
#   endif
    int Append( DataSet* );
    // ----- DataSet_Vector functions ------------ 
    void reset();
    void Resize(size_t s) {
      internalResize( s, Vec3(0.0) );
      origins_.resize( s, Vec3(0.0) );
    }
    void Resize(size_t s, Vec3 const& v) {
      internalResize(s, v);
      origins_.resize( s, Vec3(0.0) );
    }
    // -------------------------------------------

    /// \return vector origin at specified position
    const Vec3& OXYZ(int i)       const { return origins_[i]; }
    /// Add vector and origin
    void AddVxyzo(Vec3 const& v, Vec3 const& c) {
      internalAdd( v );
      origins_.push_back( c );
    }
  private:
    Varray origins_;
};
// ---------- INLINE FUNCTIONS -------------------------------------------------
void DataSet_Vector_OXYZ::Add(size_t frame, const void* vIn) {
  if (frame > Size()) {
    internalResize(frame, Vec3(0.0));
    origins_.resize(frame, Vec3(0.0));
  }
  internalAdd(vIn); 
  origins_.push_back( Vec3( (const double*)vIn+3 ) );
}
#endif
