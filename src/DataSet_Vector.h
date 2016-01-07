#ifndef INC_DATASET_VECTOR_H
#define INC_DATASET_VECTOR_H
#include "DataSet.h"
#include "DataSet_1D.h" // FIXME remove after Correlation gone
#include "Vec3.h"
#include "ComplexArray.h"
class DataSet_Vector : public DataSet {
    static const ComplexArray COMPLEXBLANK;
  public:
    static const Vec3 ZERO; ///< Vector of {0,0,0}
    typedef std::vector<Vec3> Varray;
    DataSet_Vector();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Vector();}
    // ----- DataSet functions -------------------
    size_t Size()                       const { return vectors_.size(); }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
#   endif
    void Info()                         const { return;                 }
    int Allocate(SizeArray const&);
    inline void Add(size_t, const void*);
    void WriteBuffer(CpptrajFile&,SizeArray const&) const;
    int Append( DataSet* );
    // -------------------------------------------
    void reset();
    void Resize(size_t s)               { vectors_.resize( s );    }
    void Resize(size_t s, Vec3 const& v){ vectors_.resize(s, v);   }
    //typedef Varray::iterator iterator;
    //iterator begin()                    { return vectors_.begin(); }
    //iterator end()                      { return vectors_.end();   }
    const Vec3& operator[](int i) const { return vectors_[i];      }
    Vec3&       operator[](int i)       { return vectors_[i];      }
    const Vec3& VXYZ(int i)       const { return vectors_[i];      }
    const Vec3& OXYZ(int i)       const {
      if (origins_.empty())
        return ZERO;
      return origins_[i];
    }
    bool HasOrigins()             const { return !origins_.empty(); }
    void ReserveVecs(size_t n)          { vectors_.reserve( n );   }
    void AddVxyz(Vec3 const& v)         { vectors_.push_back( v ); }
    void AddVxyz(Vec3 const& v, Vec3 const& c) {
      vectors_.push_back( v );
      origins_.push_back( c );
    }
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
  private:
    int order_;      ///< Order for spherical harmonics calculations
    Varray vectors_;
    Varray origins_;
    /// Hold spherical harmonic values for m=-order to order
    std::vector<ComplexArray> sphericalHarmonics_; // TODO Make AdditionalData
};
// ---------- INLINE FUNCTIONS -------------------------------------------------
void DataSet_Vector::Add(size_t frame, const void* vIn) {
  if (frame > vectors_.size()) {
    vectors_.resize(frame, ZERO);
    origins_.resize(frame, ZERO);
  }
  vectors_.push_back( Vec3( (const double*)vIn   ) );
  origins_.push_back( Vec3( (const double*)vIn+3 ) );
}
#endif
