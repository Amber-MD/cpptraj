#ifndef INC_DATASET_MATRIX_H
#define INC_DATASET_MATRIX_H
#include "Atom.h"
#include "AtomMask.h"
#include "DataSet.h"
/// Matrix, hold average over frames
class DataSet_Matrix : public DataSet {
  public:
    enum Matrix_Type {
      MATRIX_NULL=0, MATRIX_DIST,      MATRIX_COVAR, MATRIX_MWCOVAR,
      MATRIX_CORREL, MATRIX_DISTCOVAR, MATRIX_IDEA,  MATRIX_IRED
    };

    DataSet_Matrix();
    ~DataSet_Matrix();
    // DataSet functions -------------------------
    int Xmax() { return nrows_ - 1; }
    int Size() { return matsize_;   }
    void Write2D( CpptrajFile&, int, int);
    void GetDimensions(std::vector<int>&);
    // -------------------------------------------

    void SetType( Matrix_Type typeIn ) { type_ = typeIn; }
    void IncrementSnap() { ++snap_; }
    int Nrows() { return nrows_; }
    int Ncols() { return ncols_; }
    int Nsnap() { return snap_;  }

    int MatrixAlloc(AtomMask&, AtomMask&, std::vector<Atom> const&);
    void DivideBy(double);
    void Vect2MinusVect();

    // Iterator over matrix elements
    class iterator : public std::iterator<std::forward_iterator_tag, double>
    {
      public:
        iterator() : ptr_(0) {}
        iterator(const iterator& rhs) : ptr_(rhs.ptr_) {}
        iterator(double* pin) : ptr_(pin) {}
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
          ++ptr_;
          return *this;
        }
        iterator operator++(int) {
          iterator tmp(*this);
          ++(*this);
          return tmp;
        }
        // Value
        double& operator*() { return *ptr_; }
        // Address
        double* operator->() { return ptr_; }
        // Addition
        iterator& operator+=(int offset) {
          ptr_ += offset;
          return *this;
        }
        iterator operator+(int offset) {
          iterator tmp(*this);
          tmp += offset;
          return tmp;
        }
      private:
        double* ptr_;
    };
    iterator begin() { return mat_; }
    iterator end() { return (mat_ + matsize_); }
    iterator v1begin() { return vect_;}
    iterator v1end()   { return (vect_ + vectsize_);}
    iterator v2begin() { return vect2_;}
    iterator v2end()   { return (vect2_ + vectsize_);}
  private:
    double* mat_;     ///< Hold matrix elements
    double* vect_;    ///< Hold diagonal elements | avg coords
    double* vect2_;   ///< Square of vect_. May not need to be stored.
    double* mass_;    ///< Hold masses. Currently only for square (i.e. size is nrows)?
    int matsize_;     ///< Total number of matrix elements
    int nrows_;       ///< Number of rows in the matrix
    int ncols_;       ///< Number of columns in the matrix
    int vectsize_;    ///< Sizes of vect | vect2
    int snap_;        ///< Number of snapshots
    Matrix_Type type_; ///< Type of matrix.
    /// Pointer to index calculator for current matrix type
    int (*calcIndex)(int,int,int);

    static int calcFullIndex(int,int,int);
    static int calcHalfIndex(int,int,int);
};
#endif
