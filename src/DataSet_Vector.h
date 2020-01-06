#ifndef INC_DATASET_VECTOR_H
#define INC_DATASET_VECTOR_H
#include "DataSet.h"
#include "DataSet_1D.h" // FIXME remove after Correlation gone
#include "Vec3.h"
#include "ComplexArray.h"
class DataSet_Vector : public DataSet {
    static const ComplexArray COMPLEXBLANK;
  public:
    typedef std::vector<Vec3> Varray;
    DataSet_Vector();
    DataSet_Vector(DataType, TextFormat const&);
    // ----- DataSet functions -------------------
    size_t Size()                       const { return vectors_.size(); }
    void Info()                         const { return;                 }
    // -------------------------------------------
    virtual void reset() = 0;
    virtual void Resize(size_t) = 0;
    virtual void Resize(size_t, Vec3 const&) = 0;

    const Vec3& operator[](int i) const { return vectors_[i];      }
    Vec3&       operator[](int i)       { return vectors_[i];      }
    const Vec3& VXYZ(int i)       const { return vectors_[i];      }
    void AddVxyz(Vec3 const& v)         { vectors_.push_back( v ); }
    typedef Varray::const_iterator const_iterator;
    const_iterator begin()       const { return vectors_.begin(); }
    const_iterator end()         const { return vectors_.end();   }
    const Vec3&    Back()        const { return vectors_.back();  }

    /// Calculate auto/cross-correlation //TODO Move to Corr.cpp
    int CalcVectorCorr(DataSet_Vector const&, DataSet_1D&, int) const;
    /// Calculate spherical harmonics arrays for given Legendre order
    int CalcSphericalHarmonics(int);
    /// \return Spherical harmonics array for given m (-order_ <= m <= order_)
    ComplexArray const& SphericalHarmonics(int) const;
    /// \return Constant for normalization via spherical harmonics addition theorem.
    static double SphericalHarmonicsNorm(int);
  protected:
    size_t internalSize() const;
    void internalAlloc(SizeArray const&);
    void internalMemalloc(SizeArray const&);
    void internalCopyBlock(size_t, DataSet_Vector const&, size_t, size_t);
    void internalAppend(DataSet_Vector const*);
    void internalReset();
    inline void internalResize(size_t, Vec3 const&);
    inline void internalAdd(const void*);
    Varray const& vectors() const { return vectors_; }
#   ifdef MPI
    Varray& internalVecArray() { return vectors_; }
#   endif
  private:
    int order_;      ///< Order for spherical harmonics calculations
    Varray vectors_;
    /// Hold spherical harmonic values for m=-order to order
    std::vector<ComplexArray> sphericalHarmonics_; // TODO Make AdditionalData
};
// ---------- INLINE FUNCTIONS -------------------------------------------------
void DataSet_Vector::internalResize(size_t frame, Vec3 const& vIn) {
    vectors_.resize(frame, vIn);
}

void DataSet_Vector::internalAdd(const void* vIn) {
  vectors_.push_back( Vec3( (const double*)vIn   ) );
}
#endif
