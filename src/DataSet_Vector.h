#ifndef INC_DATASET_VECTOR_H
#define INC_DATASET_VECTOR_H
#include "DataSet_1D.h"
#include "ArgList.h"
#include "Vec3.h"
#include "ComplexArray.h"
class DataSet_Vector : public DataSet_1D {
  // TODO: Should this just be two arrays of Vec3?
  public:
    DataSet_Vector();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Vector();}
    ~DataSet_Vector();
    // ----- DataSet functions -------------------
    size_t Size()                       const { return currentidx_ / 6; }
    int Sync()                                { return 1;               }
    void Info()                         const { return; }
    // ----- DataSet_1D functions ----------------
    int Allocate1D(size_t);
    inline void Add(size_t, const void*);
    double Dval(size_t)                 const { return 0.0;             }
    void WriteBuffer(CpptrajFile&, size_t) const;
    // -------------------------------------------
    double SphereHarm(int i) const { return sphereharm_[i]; } // TODO: Replace
    void reset()                   { currentidx_ = 0;       }
    void SetIred()                 { isIred_ = true;        }
    bool IsIred()            const { return isIred_;        }
    void CalcSphericalHarmonics(int);
    int FillData(ComplexArray&, int) const; 
    // -------------------------------------------
    inline Vec3 operator[](int) const;
    inline void AddVxyz(Vec3 const&, Vec3 const&);
    inline void AddVxyz(Vec3 const&);
    inline Vec3 VXYZ(int) const;
    inline Vec3 CXYZ(int) const;
    inline Vec3 CurrentVec() const;
    // Currently only used for matrix IRED
    inline double Dot(const DataSet_Vector&) const;
    // -------------------------------------------
    // Iterator over vectors
    class iterator : public std::iterator<std::forward_iterator_tag, const double*>
    {
      public:
        iterator() : ptr_(0) {}
        iterator(const iterator& rhs) : ptr_(rhs.ptr_) {}
        iterator(const double* pin) : ptr_(pin) {}
        iterator& operator=(const iterator& rhs) {
          if (this == &rhs) return *this;
          ptr_ = rhs.ptr_;
          return *this;
        }
        // Relations
        bool operator==(const iterator& rhs) { return (ptr_==rhs.ptr_);}
        bool operator!=(const iterator& rhs) { return (ptr_!=rhs.ptr_);}
        // Increment
        iterator& operator++() {
          ptr_ += 6;
          return *this;
        }
        iterator operator++(int) {
          iterator tmp(*this);
          ++(*this);
          return tmp;
        }
        // Value
        const double* operator*() { return ptr_; }
        // Address
        const double** operator->() { return &ptr_; }
      private:
        const double* ptr_;
    };
    iterator begin() { return xyz_;               }
    iterator end()   { return xyz_ + currentidx_; }
  private:
    size_t totalidx_;    ///< Size of the xyz array
    size_t currentidx_;  ///< Current position in the xyz array
    int order_;       ///< Order of spherical harmonics
    double* xyz_;     ///< 3x Vector lengths followed by 3x vector origin
    bool writeSum_;   ///< If true will print vx+cx vy+cy vz+cz in WriteBuffer
    bool isIred_;     ///< If true this vector can be used to calc subsequent IRED matrix
    /// Array of spherical harmonic values of order order_.
    std::vector<double> sphereharm_;

    void IncreaseSize();
    static void sphericalHarmonics(int,int,const double*,double, double[2]);
};
// ---------- INLINE FUNCTIONS -------------------------------------------------
// DataSet_Vector::Add()
void DataSet_Vector::Add(size_t frameNum, const void* vIn) { 
  // TODO: Somehow incorporate frameNum?
  const double* xyz = (const double*)vIn;
  if (currentidx_ == totalidx_)
    IncreaseSize();
  xyz_[currentidx_  ] = xyz[0];
  xyz_[currentidx_+1] = xyz[1];
  xyz_[currentidx_+2] = xyz[2];
  xyz_[currentidx_+3] = xyz[3];
  xyz_[currentidx_+4] = xyz[4];
  xyz_[currentidx_+5] = xyz[5];
  currentidx_ += 6;
}
// DataSet_Vector::operator[]()
Vec3 DataSet_Vector::operator[](int i) const {
  int idx = i * 6;
  return Vec3( xyz_[idx], xyz_[idx+1], xyz_[idx+2] );
}
// DataSet_Vector::AddVxyz()
void DataSet_Vector::AddVxyz(Vec3 const& vxyz, Vec3 const& cxyz) {
  if (currentidx_ == totalidx_)
    IncreaseSize();
  xyz_[currentidx_  ] = vxyz[0];
  xyz_[currentidx_+1] = vxyz[1];
  xyz_[currentidx_+2] = vxyz[2];
  xyz_[currentidx_+3] = cxyz[0];
  xyz_[currentidx_+4] = cxyz[1];
  xyz_[currentidx_+5] = cxyz[2];
  currentidx_ += 6;
}
// DataSet_Vector::AddVxyz()
void DataSet_Vector::AddVxyz(Vec3 const& xyz) {
  if (currentidx_ == totalidx_)
    IncreaseSize();
  xyz_[currentidx_  ] = xyz[0];
  xyz_[currentidx_+1] = xyz[1];
  xyz_[currentidx_+2] = xyz[2];
  currentidx_ += 6;
}
// DataSet_Vector::VXYZ()
Vec3 DataSet_Vector::VXYZ(int i) const {
  int idx = i * 6;
  return Vec3( xyz_[idx], xyz_[idx+1], xyz_[idx+2] ); 
}
// DataSet_Vector::CXYZ()
Vec3 DataSet_Vector::CXYZ(int i) const {
  int idx = (i * 6) + 3;
  return Vec3( xyz_[idx], xyz_[idx+1], xyz_[idx+2] ); 
}
// DataSet_Vector::CurrentVec()
Vec3 DataSet_Vector::CurrentVec() const {
  if (currentidx_ == 0) return Vec3(0,0,0);
  return Vec3(xyz_[currentidx_-6], xyz_[currentidx_-5], xyz_[currentidx_-4]);
}
// DataSet_Vector::Dot()
double DataSet_Vector::Dot(const DataSet_Vector& rhs) const {
  return (xyz_[currentidx_-6]*rhs.xyz_[currentidx_-6] + 
          xyz_[currentidx_-5]*rhs.xyz_[currentidx_-5] + 
          xyz_[currentidx_-4]*rhs.xyz_[currentidx_-4]  );
}
#endif
