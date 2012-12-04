#ifndef INC_DATASET_VECTOR_H
#define INC_DATASET_VECTOR_H
#include "DataSet.h"
#include "ArgList.h"
#include "Vec3.h"
#include "ComplexArray.h"
class DataSet_Vector : public DataSet {
  // TODO: Should this just be two arrays of Vec3?
  public:
    DataSet_Vector();
    ~DataSet_Vector();
    // DataSet functions -------------------------
    int Size()  { return currentidx_ / 6;       }
    int Xmax()  { return (currentidx_ / 6) - 1; }
    int Width() { return ((width_ + 1) * 9);    }
    int Allocate(int);
    void WriteBuffer(CpptrajFile&, int);
    void Add(int frameNum, void* vIn) { // TODO: Somehow incorporate frameNum?
      double* xyz = (double*)vIn;
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
    // -------------------------------------------
    double SphereHarm(int i) { return sphereharm_[i]; } // TODO: Replace
    void reset()             { currentidx_ = 0;       }
    void SetIred()           { isIred_ = true;        }
    bool IsIred()            { return isIred_;        }
    void CalcSphericalHarmonics(int);
    int FillData(ComplexArray&, int); 
    // -------------------------------------------
    Vec3 operator[](int i) const {
      int idx = i * 6;
      return Vec3( xyz_[idx], xyz_[idx+1], xyz_[idx+2] );
    }
    void AddVxyz(Vec3 const& vxyz, Vec3 const& cxyz) {
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
    void AddVxyz(Vec3 const& xyz) {
      if (currentidx_ == totalidx_)
        IncreaseSize();
      xyz_[currentidx_  ] = xyz[0];
      xyz_[currentidx_+1] = xyz[1];
      xyz_[currentidx_+2] = xyz[2];
      currentidx_ += 6;
    }
    Vec3 VXYZ(int i) {
      int idx = i * 6;
      return Vec3( xyz_[idx], xyz_[idx+1], xyz_[idx+2] ); 
    }
    Vec3 CXYZ(int i) {
      int idx = (i * 6) + 3;
      return Vec3( xyz_[idx], xyz_[idx+1], xyz_[idx+2] ); 
    }
    Vec3 CurrentVec() {
      if (currentidx_ == 0) return Vec3(0,0,0);
      return Vec3(xyz_[currentidx_-6], xyz_[currentidx_-5], xyz_[currentidx_-4]);
    }
    // Currently only used for matrix IRED
    double Dot(const DataSet_Vector& rhs) {
      return (xyz_[currentidx_-6]*rhs.xyz_[currentidx_-6] + 
              xyz_[currentidx_-5]*rhs.xyz_[currentidx_-5] + 
              xyz_[currentidx_-4]*rhs.xyz_[currentidx_-4]  );
    }
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
    int totalidx_;    ///< Size of the xyz array
    int currentidx_;  ///< Current position in the xyz array
    int order_;       ///< Order of spherical harmonics
    std::vector<double> sphereharm_;
    double* xyz_;     ///< 3x Vector lengths followed by 3x vector origin
    bool writeSum_;   ///< If true will print vx+cx vy+cy vz+cz in WriteBuffer
    bool isIred_;     ///< If true this vector can be used to calc subsequent IRED matrix

    void IncreaseSize();
    static void sphericalHarmonics(int,int,const double*,double, double[2]);
};
#endif
