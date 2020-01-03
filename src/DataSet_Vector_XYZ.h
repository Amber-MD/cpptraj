#ifndef INC_DATASET_VECTOR_XYZ_H
#define INC_DATASET_VECTOR_XYZ_H
#include "DataSet_Vector.h"
/// Vector data set, no origins
class DataSet_Vector_XYZ : public DataSet_Vector {
  public:
    DataSet_Vector_XYZ();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Vector_XYZ();}
    // ----- DataSet functions -------------------
    size_t MemUsageInBytes() const;
    int Allocate(SizeArray const&);
    void WriteBuffer(CpptrajFile&,SizeArray const&) const;
    int Append(DataSet*);
    inline void Add(size_t, const void*);
    // ----- DataSet_Vector functions ------------
    void reset();
    void Resize(size_t s) { internalResize(s, Vec3(0.0)); }
    void Resize(size_t s, Vec3 const& v) { internalResize(s, v); }
};

// ---------- INLINE FUNCTIONS -------------------------------------------------
void DataSet_Vector_XYZ::Add(size_t frame, const void* vIn) {
  if (frame > Size())
    internalResize(frame, Vec3(0.0));
  internalAdd(vIn);
}
#endif
