#ifndef INC_DATASET_MATRIX_H
#define INC_DATASET_MATRIX_H
#include "Atom.h"
#include "AtomMask.h"
#include "DataSet.h"
#include "ArgList.h"
/// Matrix, hold average over frames
/** 
  * At the end of the calc vect will hold:
  *  -Average coordinates in the case of COVAR, MWCOVAR, CORREL;
  *  -average distances in the case of DISTCOVAR;
  *  -average of r*r/3 with r = "center of mass to atom" vector in the case of IDEA;
  *  -average of Legendre polynomial P(cos(angle(ri, rj))) in the case of IRED;
  *  -nothing in the case of DIST
  */
class DataSet_Matrix : public DataSet {
  public:
    DataSet_Matrix();
    ~DataSet_Matrix();
    // Matrix Types ------------------------------
    // NOTE: Stored here instead of in Action_Matrix since certain types
    //       of analysis depend on the type of matrix, e.g. COVAR/MWCOVAR
    //       matrices have 3 coordinates per entry as opposed to DIST which
    //       only has 1.
    enum MatrixType {
      NO_OP=0, DIST, COVAR, MWCOVAR, CORREL, DISTCOVAR, IDEA, IRED
    };
    static const char MatrixTypeString[][27];
    static const char MatrixOutputString[][10];
    static MatrixType TypeFromArg(ArgList&);
    void SetType(MatrixType typeIn) { type_ = typeIn; }
    MatrixType Type() { return type_; }
    // DataSet functions -------------------------
    int Xmax() { return nrows_ - 1; }
    int Size() { return matsize_;   }
    void Write2D( CpptrajFile&, int, int);
    void GetDimensions(std::vector<int>&);
    // -------------------------------------------

    void IncrementSnap() { ++snap_; }
    int Nrows() { return nrows_; }
    int Ncols() { return ncols_; }
    int Nelts() { return nelts_; }
    int Nsnap() { return snap_;  }
    int VectSize() { return vectsize_; }
    const double* Mass() { return mass_; }
    const double* Vect() { return vect_; }
    double* MatrixPtr()  { return mat_;  } // For interfacing w/fortran routines
    double GetElement(int col, int row) { return mat_[calcIndex(ncols_, col, row)]; } 

    int AllocateVectors(int);
    int AllocateMatrix(int,int,int,int);
    void StoreMass(std::vector<double> const&);
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
    MatrixType type_; ///< Type of matrix
    int matsize_;     ///< Total number of matrix elements
    int nrows_;       ///< Number of rows in the matrix
    int ncols_;       ///< Number of columns in the matrix
    int nelts_;       ///< Number of elements in a row of the matrix.
    int vectsize_;    ///< Sizes of vect | vect2
    int snap_;        ///< Number of snapshots
    /// Pointer to index calculator for current matrix type
    int (*calcIndex)(int,int,int);

    static int calcFullIndex(int,int,int);
    static int calcHalfIndex(int,int,int);
};
#endif
