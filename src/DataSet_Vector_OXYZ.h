#ifndef INC_DATASET_VECTOR_OXYZ_H
#define INC_DATASET_VECTOR_OXYZ_H
#include "DataSet_Vector.h"
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
      vectors_.push_back( v );
      origins_.push_back( c );
    }

    /// Calculate auto/cross-correlation //TODO Move to Corr.cpp
    int CalcVectorCorr(DataSet_Vector const&, DataSet_1D&, int) const;
    /// Calculate spherical harmonics arrays for given Legendre order
    int CalcSphericalHarmonics(int);
    /// \return Spherical harmonics array for given m (-order_ <= m <= order_)
    ComplexArray const& SphericalHarmonics(int) const;
    /// \return Constant for normalization via spherical harmonics addition theorem.
    static double SphericalHarmonicsNorm(int); 
  private:
    int order_;      ///< Order for spherical harmonics calculations
    Varray vectors_;
    Varray origins_;
    /// Hold spherical harmonic values for m=-order to order
    std::vector<ComplexArray> sphericalHarmonics_; // TODO Make AdditionalData
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
